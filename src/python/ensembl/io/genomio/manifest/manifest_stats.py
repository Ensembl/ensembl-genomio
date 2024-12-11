# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Register the main stats for all files from a manifest for comparison."""

__all__ = ["InvalidIntegrityError", "ManifestStats", "StatsLengths"]

import logging
import json
from math import floor
from os import PathLike
from pathlib import Path
from typing import Any, TypeAlias

from BCBio import GFF
from Bio import SeqIO

from ensembl.io.genomio.gff3.features import GFFSeqFeature
from ensembl.io.genomio.manifest.manifest import Manifest
from ensembl.io.genomio.utils import get_json
from ensembl.utils.archive import open_gz_file
from ensembl.utils import StrPath


# Record the lengths of the sequence for features/regions
StatsLengths: TypeAlias = dict[str, int]


class InvalidIntegrityError(Exception):
    """When a file integrity check fails"""


class ManifestStats:
    """Representation of the main stats of the files in a manifest for comparison.

    The stats in question are:
    - lengths of sequences (DNA, genes and peptides)
    - sequences and features IDs
    - sequences circularity
    """

    def __init__(self, manifest_path: StrPath, ignore_final_stops: bool = False) -> None:
        self.manifest_files = self._get_manifest(manifest_path)
        self.genome: dict[str, Any] = {}

        self.lengths: dict[str, StatsLengths] = {
            "dna_sequences": {},
            "peptide_sequences": {},
            "seq_region_levels": {},
            "annotations": {},
            "agp": {},
            "gff3_seq_regions": {},
            "gff3_genes": {},
            "gff3_translations": {},
            "gff3_all_translations": {},
            "gff3_transposable_elements": {},
            "ann_genes": {},
            "ann_translations": {},
            "ann_transposable_elements": {},
            "seq_regions": {},
        }

        self.circular: dict[str, StatsLengths] = {
            "seq_regions": {},
        }

        self.errors: list[str] = []

        self.ignore_final_stops = ignore_final_stops

    def _get_manifest(self, manifest_path: PathLike) -> dict[str, Any]:
        """Load the content of a manifest file.

        Returns:
            Dict: Content of the manifest file.
        """
        manifest = Manifest(Path(manifest_path).parent)
        manifest_files = manifest.load()

        # Replace the {file, md5} dict with the file path
        for name in manifest_files:
            if "file" in manifest_files[name]:
                manifest_files[name] = Path(manifest_path).parent / manifest_files[name]["file"]
            else:
                for f in manifest_files[name]:
                    manifest_files[name][f] = Path(manifest_path).parent / manifest_files[name][f]["file"]
        return manifest_files

    def add_error(self, error: str) -> None:
        """Record an error."""
        self.errors.append(error)

    def load_seq_regions(self) -> None:
        """Retrieve seq_regions lengths and circular information from the seq_region JSON file."""

        if "seq_region" not in self.manifest_files:
            return
        logging.info("Manifest contains seq_region JSON")
        seq_regions = get_json(Path(self.manifest_files["seq_region"]))
        seqr_seqlevel = {}
        seq_lengths = {}
        seq_circular = {}
        # Store the length as int
        for seq in seq_regions:
            seq_lengths[seq["name"]] = int(seq["length"])
            seq_circular[seq["name"]] = seq.get("circular", False)
            if seq["coord_system_level"] == "contig":
                seqr_seqlevel[seq["name"]] = int(seq["length"])
            # Also record synonyms (in case GFF file uses synonyms)
            if "synonyms" in seq:
                for synonym in seq["synonyms"]:
                    seq_lengths[synonym["name"]] = int(seq["length"])
        self.lengths["seq_regions"] = seq_lengths
        self.circular["seq_regions"] = seq_circular

    def load_peptides_fasta_lengths(self) -> None:
        """Retrieve peptides sequences lengths from their FASTA file."""
        if "fasta_pep" not in self.manifest_files:
            return
        self.lengths["peptide_sequences"] = self._get_fasta_lengths(
            self.manifest_files["fasta_pep"], ignore_final_stops=self.ignore_final_stops
        )

    def load_dna_fasta_lengths(self) -> None:
        """Retrieve DNA sequences lengths from their FASTA file."""
        if "fasta_dna" not in self.manifest_files:
            return
        self.lengths["dna_sequences"] = self._get_fasta_lengths(self.manifest_files["fasta_dna"])

    def _get_fasta_lengths(self, fasta_path: StrPath, ignore_final_stops: bool = False) -> dict[str, int]:
        """Returns every sequence ID and its length from a FASTA file (DNA or peptide).

        An error will be added for every empty id, non-unique id or stop codon found in the FASTA file.

        Args:
            fasta_path: Path to FASTA file.
            ignore_final_stops: Do not include final stop in the total length.

        """

        data = {}
        non_unique = {}
        non_unique_count = 0
        empty_id_count = 0
        contains_stop_codon = 0
        rec_count = 0
        for rec in SeqIO.parse(fasta_path, "fasta"):
            rec_count += 1

            # Flag empty ids
            if rec.id == "":
                empty_id_count += 1
                continue
            # Flag redundant ids
            if rec.id in data:
                non_unique[rec.id] = 1
                non_unique_count += 1
            # Store sequence id and length
            data[rec.id] = len(rec.seq)
            stops = rec.seq.count("*")
            if stops >= 1 and not rec.seq.endswith("*"):
                contains_stop_codon += 1
            elif rec.seq.endswith("*") and not ignore_final_stops:
                contains_stop_codon += 1

        if empty_id_count > 0:
            self.add_error(f"{empty_id_count} sequences with empty ids in {fasta_path}")
        if non_unique_count > 0:
            self.add_error(f"{non_unique_count} non unique sequence ids in {fasta_path}")
        if contains_stop_codon > 0:
            self.add_error(f"{contains_stop_codon} sequences with stop codons in {fasta_path}")
        if rec_count == 0:
            self.add_error(f"No sequences found in {fasta_path}")
        return data

    def load_functional_annotation(self) -> None:
        """Load the functional annotation file to retrieve the gene_id and translation id.

        The functional annotation file is stored in a JSON format containing the description, id
        and object type, eg: "gene", "transcript", "translation".

        """
        if "functional_annotation" not in self.manifest_files:
            return
        logging.info("Manifest contains functional annotation(s)")

        # Load the json file
        with open(self.manifest_files["functional_annotation"]) as json_file:
            data = json.load(json_file)

        # Get gene ids and translation ids
        genes = {}
        translations = {}
        transposons = {}

        for item in data:
            if item["object_type"] == "gene":
                genes[item["id"]] = 1
            elif item["object_type"] == "translation":
                translations[item["id"]] = 1
            if item["object_type"] == "transposable_element":
                transposons[item["id"]] = 1

        stats = {
            "ann_genes": genes,
            "ann_translations": translations,
            "ann_transposable_elements": transposons,
        }
        self.lengths = {**self.lengths, **stats}

    def load_gff3(self) -> None:
        """A GFF3 parser is used to retrieve information in the GFF3 file such as
        gene and CDS ids and their corresponding lengths.
        """
        if "gff3" not in self.manifest_files:
            return
        logging.info("Manifest contains GFF3 gene annotations")
        gff3_path = self.manifest_files["gff3"]

        seqs: StatsLengths = {}
        genes: StatsLengths = {}
        peps: StatsLengths = {}
        all_peps: StatsLengths = {}
        tes: StatsLengths = {}

        with open_gz_file(gff3_path) as gff3_handle:
            gff = GFF.parse(gff3_handle)
            for seq in gff:
                seqs[seq.id] = len(seq.seq)

                for feat in seq.features:
                    feat_length = abs(feat.location.end - feat.location.start)
                    # Store gene id and length
                    if feat.type in ["gene", "ncRNA_gene", "pseudogene"]:
                        self._retrieve_gff_gene_lengths(feat, genes, peps, all_peps)
                    if feat.type == "transposable_element":
                        tes[feat.id] = feat_length

        stats: dict[str, StatsLengths] = {
            "gff3_seq_regions": seqs,
            "gff3_genes": genes,
            "gff3_translations": peps,
            "gff3_all_translations": all_peps,
            "gff3_transposable_elements": tes,
        }
        self.lengths = {**self.lengths, **stats}

    def _retrieve_gff_gene_lengths(
        self, feat: GFFSeqFeature, genes: StatsLengths, peps: StatsLengths, all_peps: StatsLengths
    ) -> None:
        """Record genes and peptides lengths from a feature.

        Args:
            feat : Gene feature to check.
            genes: Record of genes lengths to update.
            peps: Record of peptides lengths to update.
            all_peps: Record of all peptides lengths to update (include pseudogenes).

        """
        gene_id = feat.id
        gene_id = gene_id.replace("gene:", "")
        genes[gene_id] = abs(feat.location.end - feat.location.start)
        # Get CDS id and length
        protein_transcripts = {
            "mRNA",
            "pseudogenic_transcript",
        }
        ig_transcripts = {
            "IG_V_gene",
            "IG_C_gene",
            "TR_C_gene",
            "TR_V_gene",
        }
        cds_transcripts = protein_transcripts.union(ig_transcripts)
        for feat2 in feat.sub_features:
            if feat2.type not in cds_transcripts:
                continue
            length: dict[str, int] = {}
            for feat3 in feat2.sub_features:
                if feat3.type != "CDS":
                    continue
                pep_id = feat3.id
                pep_id = pep_id.replace("CDS:", "")
                length.setdefault(pep_id, 0)
                length[pep_id] += abs(feat3.location.end - feat3.location.start)
            for pep_id, pep_length in length.items():
                # Store length for translations, add pseudo translations separately
                pep_length = floor(pep_length / 3) - 1
                if feat.type != "pseudogene" and feat2.type in protein_transcripts:
                    peps[pep_id] = pep_length
                all_peps[pep_id] = pep_length

    def load_agp_seq_regions(self, agp_dict: dict | None) -> None:
        """AGP files describe the assembly of larger sequence objects using smaller objects.

        E.g. describes the assembly of scaffolds from contigs.

        Args:
            agp_dict: Dict containing the information about the sequence.

        Note:
            AGP file is only used in the older builds, not used for current processing.
        """
        if not agp_dict:
            return
        logging.info("Manifest contains AGP files")

        seqr: StatsLengths = {}
        for agp_path in agp_dict.values():
            with open(agp_path, "r") as agph:
                for line in agph:
                    (
                        asm_id,
                        _,  # asm_start
                        asm_end,
                        _,  # asm_part
                        typ,
                        cmp_id,
                        _,  # cmp_start
                        cmp_end,
                        _,  # cmp_strand
                    ) = line.split("\t")
                    # Ignore WGS contig
                    if typ != "W":
                        continue

                    # Assembled seq length
                    if asm_id not in seqr or seqr[asm_id] < int(asm_end):
                        seqr[asm_id] = int(asm_end)

                    # Composite seq length
                    if cmp_id not in seqr or seqr[cmp_id] < int(cmp_end):
                        seqr[cmp_id] = int(cmp_end)

        self.lengths["agp"] = seqr

    def load_genome(self) -> None:
        """Load the genome data."""
        if "genome" not in self.manifest_files:
            return
        logging.info("Manifest contains genome JSON")
        self.genome = get_json(Path(self.manifest_files["genome"]))

    def prepare_integrity_data(self) -> None:  # pylint: disable=too-many-branches
        """Read all the files and keep a record (IDs and their lengths) for each case to be compared later."""
        self.load_gff3()
        self.load_dna_fasta_lengths()
        self.load_peptides_fasta_lengths()
        self.load_seq_regions()
        self.load_functional_annotation()
        self.load_agp_seq_regions(self.manifest_files.get("agp"))
        self.load_genome()

    def has_lengths(self, name: str) -> bool:
        """Check if a given name has lengths records.

        Args:
            name: Name for the lengths to check.

        Raises:
            KeyError: If the name is not supported.
        """
        try:
            return bool(self.lengths[name])
        except KeyError as err:
            raise KeyError(f"There is no length record for {name}") from err

    def get_lengths(self, name: str) -> dict[str, Any]:
        """Returns a dict associating IDs with their length from a given file name.

        Args:
            name: Name for the lengths to get.

        Raises:
            KeyError: If the name is not supported.
        """
        try:
            return self.lengths[name]
        except KeyError as err:
            raise KeyError(f"There is no length record for {name}") from err

    def get_circular(self, name: str) -> dict[str, Any]:
        """Returns a dict associating IDs with their is_circular flag from a given file name.

        Args:
            name: Name for the circular data to get.

        Raises:
            KeyError: If the name is not supported.
        """
        try:
            return self.circular[name]
        except KeyError as err:
            raise KeyError(f"No length available for key {name}") from err
