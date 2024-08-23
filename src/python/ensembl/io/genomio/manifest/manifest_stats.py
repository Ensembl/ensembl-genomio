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
from typing import Any

from BCBio import GFF
from Bio import SeqIO

from ensembl.io.genomio.utils import get_json
from ensembl.io.genomio.gff3.features import GFFSeqFeature
from ensembl.io.genomio.manifest.manifest import Manifest
from ensembl.utils.archive import open_gz_file


# Record the lengths of the sequence for features/regions
StatsLengths = dict[str, int]


class InvalidIntegrityError(Exception):
    """When a file integrity check fails"""


class ManifestStats:
    """Representation of the main stats of the files in a manifest for comparison.

    The stats in question are:
    - lengths of sequences (DNA, genes and peptides)
    - sequences and features IDs
    - sequences circularity
    """

    def __init__(self, manifest_path: PathLike) -> None:
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

        self.ignore_final_stops = False
        self.brc_mode = False

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

    def get_seq_regions(self):
        """Retrieve seq_regions lengths and circular information from the seq_region JSON file."""
        if not "seq_region" in self.manifest_files:
            return
        logging.info("Manifest contains seq_region JSON")
        seq_regions = get_json(Path(self.manifest_files["seq_region"]))
        if len(seq_regions) == 0:
            self.add_error(f"No sequences found in {self.manifest_files['seq_region']}")
            return
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

    def get_fasta_lengths(self, fasta_path, ignore_final_stops=False):
        """Check if the fasta files have the correct ids and no stop codon.

        Args:
            fasta_path: Path to fasta_dna and fasta_pep files.

        Returns:
            Error if any empty ids, non-unique ids or stop codons are found in the fasta files.
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
            else:
                # Flag redundant ids
                if rec.id in data:
                    non_unique[rec.id] = 1
                    non_unique_count += 1
                # Store sequence id and length
                data[rec.id] = len(rec.seq)
                stops = rec.seq.count("*")
                if stops > 1:
                    contains_stop_codon += 1
                elif stops == 1:
                    if not rec.seq.endswith("*") or not ignore_final_stops:
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

    def get_functional_annotation(self, json_path: Path | None) -> None:
        """Load the functional annotation file to retrieve the gene_id and translation id.
            A functional annotation file contains information about a gene.
            The functional annotation file is stored in a json format containing
            the description, id and object type (eg: "gene", "transcript", "translation").

        Args:
            json_path: Path to functional_annotation.json.
        """
        if not json_path:
            return
        logging.info("Manifest contains functional annotation(s)")

        # Load the json file
        with open(json_path) as json_file:
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

    def get_gff3(self, gff3_path: Path) -> None:
        """A GFF parser is used to retrieve information in the GFF file such as
           gene and CDS ids and their corresponding lengths.

        Args:
            gff3_path: Path to gff3 file.
        """

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
        if not self.brc_mode:
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
            length = {}
            for feat3 in feat2.sub_features:
                if feat3.type != "CDS":
                    continue
                pep_id = feat3.id
                if not self.brc_mode:
                    pep_id = pep_id.replace("CDS:", "")
                if pep_id not in length:
                    length[pep_id] = 0
                length[pep_id] += abs(feat3.location.end - feat3.location.start)
            for pep_id, pep_length in length.items():
                # Store length for translations, add pseudo translations separately
                pep_length = floor(pep_length / 3) - 1
                if feat.type != "pseudogene" and feat2.type in protein_transcripts:
                    peps[pep_id] = pep_length
                all_peps[pep_id] = pep_length

    def get_agp_seq_regions(self, agp_dict: dict | None):
        """AGP files describe the assembly of larger sequence objects using smaller objects.
            Eg: describes the assembly of scaffolds from contigs.

        Args:
            agp_dict: dict containing the information about the sequence.

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

    def prepare_integrity_data(self) -> None:  # pylint: disable=too-many-branches
        """Read all the files and keep a record (IDs and their lengths)
        for each cases to be compared later.
        """
        # First, get the Data
        if "gff3" in self.manifest_files:
            logging.info("Manifest contains GFF3")
            self.get_gff3(self.manifest_files["gff3"])
        if "fasta_dna" in self.manifest_files:
            logging.info("Manifest contains DNA fasta")
            # Verify if the length and id for the sequence is unique
            self.lengths["dna_sequences"] = self.get_fasta_lengths(self.manifest_files["fasta_dna"])
        if "fasta_pep" in self.manifest_files:
            logging.info("Manifest contains Peptide fasta")
            # Verify if the length and id for the sequence is unique
            self.lengths["peptide_sequences"] = self.get_fasta_lengths(
                self.manifest_files["fasta_pep"], ignore_final_stops=self.ignore_final_stops
            )
            self.get_seq_regions()
        self.get_functional_annotation(self.manifest_files.get("functional_annotation"))
        self.get_agp_seq_regions(self.manifest_files.get("agp"))
        if "genome" in self.manifest_files:
            logging.info("Manifest contains genome JSON")
            self.lengths["genome"] = get_json(Path(self.manifest_files["genome"]))

    def has_lengths(self, name: str) -> bool:
        """Check if a given name has lengths records.

        Raise KeyError if the name is not supported.

        """
        try:
            return bool(self.lengths[name])
        except KeyError as err:
            raise KeyError(f"There is no length record for {name}") from err

    def get_lengths(self, name: str) -> dict[str, Any]:
        """Returns a dict associating IDs with their length from a given file name."""
        try:
            return self.lengths[name]
        except KeyError as err:
            raise KeyError(f"There is no length record for {name}") from err

    def get_circular(self, name: str) -> dict[str, Any]:
        """Returns a dict associating IDs with their is_circular flag from a given file name."""
        try:
            return self.circular[name]
        except KeyError as err:
            raise KeyError(f"No length available for key {name}") from err
