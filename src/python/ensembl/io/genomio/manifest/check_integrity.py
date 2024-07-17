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
"""Compare the genomic data in a DNA fasta file, seq_region JSON, gene models GFF3 and peptide fasta
to ensure their contents are in sync.
"""

__all__ = ["InvalidIntegrityError", "Manifest", "IntegrityTool"]

import hashlib
import logging
import json
from math import floor
from os import PathLike
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Union

from BCBio import GFF
from Bio import SeqIO

from ensembl.io.genomio.utils import get_json
from ensembl.io.genomio.gff3.features import GFFSeqFeature
from ensembl.utils.archive import open_gz_file
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


# Record the lengths of the sequence for features/regions
Lengths = Dict[str, int]


class InvalidIntegrityError(Exception):
    """When a file integrity check fails"""


class Manifest:
    """Representation of the manifest and its files."""

    def __init__(self, manifest_path: PathLike) -> None:
        self.manifest_files = self.get_manifest(manifest_path)
        self.genome: Dict[str, Any] = {}
        self.seq_regions: Dict[str, Any] = {}

        self.lengths: Dict[str, Lengths] = {
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

        self.circular: Dict[str, Lengths] = {
            "seq_regions": {},
        }

        self.errors: List[str] = []

        self.ignore_final_stops = False
        self.brc_mode = False

    def has_lengths(self, name: str) -> bool:
        """Check if a given name has lengths records.

        Raise KeyError if the name is not supported.

        """
        try:
            if self.lengths[name]:
                return True
            return False
        except KeyError as err:
            raise KeyError(f"There is no length record for {name}") from err

    def get_lengths(self, name: str) -> Dict[str, Any]:
        """Returns a dict associating IDs with their length from a given file name."""
        try:
            return self.lengths[name]
        except KeyError as err:
            raise KeyError(f"There is no length record for {name}") from err

    def get_circular(self, name: str) -> Dict[str, Any]:
        """Returns a dict associating IDs with their is_circular flag from a given file name."""
        try:
            return self.circular[name]
        except KeyError as err:
            raise KeyError(f"No length available for key {name}") from err

    def _add_error(self, error: str) -> None:
        self.errors.append(error)

    def get_manifest(self, manifest_path: PathLike) -> Dict[str, Any]:
        """Load the content of a manifest file.

        Returns:
            Dict: Content of the manifest file.
        """
        manifest_path = Path(manifest_path)
        with manifest_path.open("r") as manifest_fh:
            manifest = json.load(manifest_fh)

            # Use dir name from the manifest
            for name in manifest:
                if "file" in manifest[name]:
                    file_path = manifest_path.parent / manifest[name]["file"]
                    # check if the md5sum is correct
                    md5sum = manifest[name]["md5sum"]
                    self._check_md5sum(file_path, md5sum)

                    manifest[name] = file_path
                else:
                    for f in manifest[name]:
                        if "file" in manifest[name][f]:
                            file_path = manifest_path.parent / manifest[name][f]["file"]
                            # check if the md5sum is correct
                            md5sum = manifest[name][f]["md5sum"]
                            self._check_md5sum(file_path, md5sum)

                            manifest[name][f] = file_path
            return manifest

    def _check_md5sum(self, file_path: Path, md5sum: str) -> None:
        """Verify the integrity of the files in manifest.json.

        An MD5 hash is generated using the path provided which is then compared to the hash in manifest.json.
        Errors are stored in ``self.errors``.

        Args:
            file_path: Path to a genome file.
            md5sum: MD5 hash for the files.
        """

        with file_path.open("rb") as f:
            bytes_obj = f.read()
            readable_hash = hashlib.md5(bytes_obj).hexdigest()
            if readable_hash != md5sum:
                raise InvalidIntegrityError(f"Invalid md5 checksum for {file_path}")

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
        if "seq_region" in self.manifest_files:
            logging.info("Manifest contains seq_region JSON")
            seq_regions = get_json(Path(self.manifest_files["seq_region"]))
            if len(seq_regions) == 0:
                self._add_error(f"No sequences found in {self.manifest_files['seq_region']}")
            else:
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
                self.seq_regions = seq_regions
        if "functional_annotation" in self.manifest_files:
            logging.info("Manifest contains functional annotation(s)")
            self.get_functional_annotation(self.manifest_files["functional_annotation"])
        if "agp" in self.manifest_files:
            logging.info("Manifest contains AGP files")
            self.lengths["agp"] = self.get_agp_seq_regions(self.manifest_files["agp"])
        if "genome" in self.manifest_files:
            logging.info("Manifest contains genome JSON")
            self.lengths["genome"] = get_json(Path(self.manifest_files["genome"]))

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
            self._add_error(f"{empty_id_count} sequences with empty ids in {fasta_path}")
        if non_unique_count > 0:
            self._add_error(f"{non_unique_count} non unique sequence ids in {fasta_path}")
        if contains_stop_codon > 0:
            self._add_error(f"{contains_stop_codon} sequences with stop codons in {fasta_path}")
        if rec_count == 0:
            self._add_error(f"No sequences found in {fasta_path}")
        return data

    def get_functional_annotation(self, json_path: Path) -> None:
        """Load the functional annotation file to retrieve the gene_id and translation id.
            A functional annotation file contains information about a gene.
            The functional annotation file is stored in a json format containing
            the description, id and object type (eg: "gene", "transcript", "translation").

        Args:
            json_path: Path to functional_annotation.json.

        Returns:
            dict with gene and translation ids.
        """

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

        seqs: Lengths = {}
        genes: Lengths = {}
        peps: Lengths = {}
        all_peps: Lengths = {}
        tes: Lengths = {}

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

        stats: Dict[str, Lengths] = {
            "gff3_seq_regions": seqs,
            "gff3_genes": genes,
            "gff3_translations": peps,
            "gff3_all_translations": all_peps,
            "gff3_transposable_elements": tes,
        }
        self.lengths = {**self.lengths, **stats}

    def _retrieve_gff_gene_lengths(
        self, feat: GFFSeqFeature, genes: Lengths, peps: Lengths, all_peps: Lengths
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
            if feat2.type in cds_transcripts:
                length = {}
                for feat3 in feat2.sub_features:
                    if feat3.type == "CDS":
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

    def get_agp_seq_regions(self, agp_dict):
        """AGP files describe the assembly of larger sequence objects using smaller objects.
            Eg: describes the assembly of scaffolds from contigs.

        Args:
            agp_dict: dict containing the information about the sequence.

        Note:
            AGP file is only used in the older builds, not used for current processing.
        """

        seqr = {}
        for agp in agp_dict:
            agp_path = agp_dict[agp]

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

        return seqr


class IntegrityTool:
    """Check the integrity of sequence and annotation files in the genome"""

    def __init__(
        self,
        manifest_file: Path,
        brc_mode: bool = False,
        ignore_final_stops: bool = False,
        no_fail: bool = False,
    ) -> None:
        self.manifest = Manifest(manifest_file)
        self.brc_mode = False
        self.set_brc_mode(brc_mode)
        self.ignore_final_stops = False
        self.set_ignore_final_stops(ignore_final_stops)
        self.errors: List[str] = []
        self.no_fail = no_fail

    def add_errors(self, errors: Union[List[str], str]) -> None:
        """Store the given errors (list or single string) in the list of all errors."""
        if isinstance(errors, str):
            self.errors.append(errors)
        else:
            self.errors += errors

    def check_integrity(self):
        """Load files listed in the manifest.json and check the integrity.
        Check if the files are correct by verifying the MD5 hash.
        Check if translation, functional annotation and sequence region ids
        and lengths are consistent with the information in gff.
        Compare sequence length from fasta_dna file to seq_region.json metadata.
        """

        # Load the manifest integrity counts
        manifest = self.manifest
        manifest.prepare_integrity_data()

        genome = manifest.genome
        seq_regions = manifest.seq_regions

        dna = manifest.get_lengths("dna_sequences")
        pep = manifest.get_lengths("peptide_sequences")
        seq_lengths = manifest.get_lengths("seq_regions")
        seq_circular = manifest.get_circular("seq_regions")

        agp_seqr = manifest.get_lengths("agp")

        # Then, run the checks
        self._check_genome(genome)

        if self.manifest.errors:
            errors_str = "\n".join(self.manifest.errors)
            raise InvalidIntegrityError(f"Manifest files parsing failed:\n{errors_str}")

        # Check gff3
        if manifest.has_lengths("gff3_genes"):
            gff_genes = manifest.get_lengths("gff3_genes")
            gff_seq_regions = manifest.get_lengths("gff3_seq_regions")
            gff_translations = manifest.get_lengths("gff3_translations")
            gff_all_translations = manifest.get_lengths("gff3_all_translations")
            gff_transposable_elements = manifest.get_lengths("gff3_transposable_elements")

            ann_genes = manifest.get_lengths("ann_genes")
            ann_translations = manifest.get_lengths("ann_translations")
            ann_transposable_elements = manifest.get_lengths("ann_transposable_elements")

            # Check fasta_pep.fa integrity
            # The sequence length and id retrieved from the fasta_pep file
            # and compared to the translated CDS id and length in the gff
            # We do not compare the peptide lengths because of sequence edits
            if pep:
                tr_errors = self.check_lengths(
                    pep, gff_translations, "Fasta translations vs gff", special_diff=True
                )
                if len(tr_errors) > 0:
                    # The pseudo CDSs are included in this check
                    # Pseudo CDSs are not translated, if the pseudo translation ids are not ignored
                    # in the gff it will give an error
                    tr_errors_all = self.check_lengths(
                        pep,
                        gff_all_translations,
                        "Fasta translations vs gff (include pseudo CDS)",
                        special_diff=True,
                    )
                    if tr_errors_all:
                        self.add_errors(tr_errors)
                        self.add_errors(tr_errors_all)

            # Check functional_annotation.json integrity
            # Gene ids, translated CDS ids and translated CDSs
            # including pseudogenes are compared to the gff
            if ann_genes:
                self.add_errors(self.check_ids(ann_genes, gff_genes, "Gene ids metadata vs gff"))
                tr_id_errors = self.check_ids(
                    ann_translations, gff_translations, "Translation ids metadata vs gff"
                )
                if tr_id_errors:
                    tr_id_errors_all = self.check_ids(
                        ann_translations,
                        gff_all_translations,
                        "Translation ids metadata vs gff (include pseudo CDS)",
                    )
                    if tr_id_errors_all:
                        self.add_errors(tr_id_errors)
                        self.add_errors(tr_id_errors_all)
                self.add_errors(
                    self.check_ids(
                        ann_transposable_elements,
                        gff_transposable_elements,
                        "TE ids metadata vs gff",
                    )
                )

            # Check the seq.json integrity
            # Compare the length and id retrieved from seq.json to the gff
            if seq_regions:
                self.check_seq_region_lengths(
                    seq_lengths, gff_seq_regions, "Seq_regions metadata vs gff", seq_circular
                )

        # Check fasta dna and seq_region integrity
        if dna and seq_regions:
            self.check_seq_region_lengths(seq_lengths, dna, "seq_regions json vs dna")

        # Check agp and seq_region integrity
        if agp_seqr and seq_lengths:
            self.check_seq_region_lengths(seq_lengths, agp_seqr, "seq_regions json vs agps")

        if self.errors:
            errors_str = "\n".join(self.errors)
            message = f"Integrity test failed:\n{errors_str}"
            if self.no_fail:
                print(message)
            else:
                raise InvalidIntegrityError(message)

    def set_brc_mode(self, brc_mode: bool) -> None:
        """Set brc mode for this tool and the manifest."""
        self.brc_mode = brc_mode
        self.manifest.brc_mode = brc_mode

    def set_ignore_final_stops(self, ignore_final_stops: bool) -> None:
        """Set ignore_final_stops (when calculating peptide length) for this tool and the manifest."""
        self.ignore_final_stops = ignore_final_stops
        self.manifest.ignore_final_stops = ignore_final_stops

    def _check_genome(self, genome: Dict) -> None:
        """Check if the accession is correct in genome.json."""
        if genome:
            if "assembly" in genome:
                genome_ass = genome["assembly"]
                if "accession" in genome_ass:
                    genome_acc = genome_ass["accession"]
                    if not re.match(r"GC[AF]_\d{9}(\.\d+)?", genome_acc):
                        self.add_errors(f"Genome assembly accession is wrong: '{genome_acc}'")

    def check_ids(self, list1, list2, name) -> List[str]:
        """Compare the ids in list1 and list2.

        Args:
            list1: dict containing sequence ids retrieved from functional.json.
            list2: dict containing length and id in the retrieved from the gff.
            name:  string

        Return:
            Whether the checks found errors.
        """

        only1 = []
        only2 = []
        common = []

        for item_id in list1:
            if item_id in list2:
                common.append(item_id)
            else:
                only1.append(item_id)
        for item_id in list2:
            if item_id not in common:
                only2.append(item_id)

        errors = []
        if common:
            logging.info(f"{len(common)} common elements in {name}")
        if only1:
            errors.append(f"{len(only1)} only in first list in {name} (first: {only1[0]})")
            logging.debug(f"{len(only1)} only in first list in {name}")
        if only2:
            errors.append(f"{len(only2)} only in second list in {name} (first: {only2[0]})")
            logging.debug(f"{len(only1)} only in second list in {name}")

        return errors

    def check_lengths(self, list1, list2, name, allowed_len_diff=None, special_diff=False) -> List[str]:
        """Check the difference in ids and length between list1 and list2.
            There are a few special cases here where we allow a certain asymmetry
            by changing the values of the arguments.

        Args:
            list1: dict containing length and id of the sequence from fasta files.
            list2: dict containing length and id in the retrieved from the gff.
            name:  string

        allowed_len_diff : None to to not accept differences in length between list1 and list2.
            The value can be changed based on how much difference in sequence length we are wanting to accept.

        special_diff: set as False when no special length difference is expected between the lists.
                    This can be changed if we want to report common sequences with 1 BP difference.

        Returns:
            Error if there is a difference in length or ids between the lists.
        """

        # check list differences, checks if abs(values diff) < allowed_len_diff

        set1 = frozenset(list1)
        set2 = frozenset(list2)
        list1_2 = list(set1 - set2)
        list2_1 = list(set2 - set1)

        errors = []
        if len(list1_2) > 0:
            errors.append(f"{name}: {len(list1_2)} from the first list only (i.e. {list1_2[0]})")
        if len(list2_1) > 0:
            errors.append(f"{name}: {len(list2_1)} from the second list only (i.e. {list2_1[0]})")

        common_len = 0
        if allowed_len_diff is None:
            common_len = len(set1 & set2)
        else:
            # check for the sequence length difference
            diff_len_list: List[str] = []
            diff_len_special_list: List[str] = []
            for e in set1 & set2:
                dl12 = list1[e] - list2[e]
                if abs(dl12) <= allowed_len_diff:
                    common_len += 1
                else:
                    _dlist = diff_len_list
                    # Special case: 1 AA /BP shorter,
                    #   so assuming the stop codon is not included in the CDS (when it should be)
                    if dl12 == 1 and special_diff:
                        _dlist = diff_len_special_list
                    _dlist.append(f"{e}: {list1[e]}, {list2[e]}")
            if diff_len_special_list:
                errors.append(
                    (
                        f"{len(diff_len_special_list)} common elements with one BP/AA length diff for {name}"
                        f"(e.g. {diff_len_special_list[0]})"
                    )
                )
            if diff_len_list:
                errors.append(
                    (
                        f"{len(diff_len_list)} common elements with length diff for {name}"
                        f"(e.g. {diff_len_list[0]})"
                    )
                )
        if common_len > 0:
            logging.warning(f"{common_len} common elements between lists for {name}")

        return errors

    def check_seq_region_lengths(
        self,
        seqrs: Dict[str, Any],
        feats: Dict[str, Any],
        name: str,
        circular: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Check the integrity of seq_region.json file by comparing the length of the sequence
            to fasta files and the gff.

            Seq_region file is in json format containing the metadata of the sequence.
            It contains sequence id, length, location and the synonyms for the sequence name
            from different sources.

        Args:
            seqs: Sequence name and length retrieved from seq_region.json file.
            feats: Sequence name and length retrieved from the fasta and gff file.
            name: Name of the check to show in the logs.
            circular: Whether any sequence is circular.

        Returns:
            Error if there are common sequences with difference in ids
            and if the sequences are not consistent in the files.
        """
        comp = self._compare_seqs(seqrs, feats, circular)

        common = comp["common"]
        diff = comp["diff"]
        diff_circular = comp["diff_circular"]
        only_seqr = comp["only_seqr"]
        only_feat = comp["only_feat"]

        if common:
            logging.info(f"{len(common)} common elements in {name}")
        if diff_circular:
            example = diff_circular[0]
            logging.info(f"{len(diff_circular)} differences for circular elements in {name} (e.g. {example})")
        if diff:
            self.add_errors(f"{len(diff)} common elements with higher length in {name} (e.g. {diff[0]})")
        if only_seqr:
            # Not an error!
            logging.info(f"{len(only_seqr)} only in seq_region list in {name} (first: {only_seqr[0]})")
        if only_feat:
            self.add_errors(f"{len(only_feat)} only in second list in {name} (first: {only_feat[0]})")

    def _compare_seqs(
        self, seqrs: Dict[str, Any], feats: Dict[str, Any], circular: Optional[Dict[str, Any]] = None
    ) -> Dict[str, List[str]]:
        """Give the intersection and other comparison between two groups of sequences.

        Args:
            seqs: Sequence name and length retrieved from seq_region.json file.
            feats: Sequence name and length retrieved from the fasta and gff file.
            circular: Whether any sequence is circular.

        Returns: Dict with 5 stats:
            common: Common elements.
            only_seqr: Elements only in the first one.
            only_feat: Elements only in the second one.
            diff: Elements that differ.
            diff_circular: Elements that differ in a circular sequence.

        """
        comp: Dict[str, List[str]] = {
            "common": [],
            "only_seqr": [],
            "only_feat": [],
            "diff": [],
            "diff_circular": [],
        }

        for seq_id in seqrs:
            if seq_id in feats:
                # Check that feature is within the seq_region length
                if feats[seq_id] > seqrs[seq_id]:
                    diff_str = f"{seq_id}: {seqrs[seq_id]} vs {feats[seq_id]}"
                    if circular and circular.get(seq_id, False):
                        comp["diff_circular"].append(diff_str)
                    else:
                        comp["diff"].append(diff_str)
                else:
                    comp["common"].append(seq_id)
            else:
                comp["only_seqr"].append(seq_id)

        for seq_id in feats:
            if (
                seq_id not in comp["common"]
                and seq_id not in comp["diff"]
                and seq_id not in comp["diff_circular"]
                and seq_id not in seqrs
            ):
                comp["only_feat"].append(seq_id)

        return comp


def main() -> None:
    """Main entrypoint."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument_src_path("--manifest_file", required=True, help="Manifest file for the data to check")
    parser.add_argument("--brc_mode", action="store_true", help="Enable BRC mode")
    parser.add_argument(
        "--ignore_final_stops", action="store_true", help="Ignore final stop when calculating peptide length"
    )
    parser.add_argument(
        "--no_fail", action="store_true", help="In case of errors, don't fail but print errors to stdout."
    )
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    inspector = IntegrityTool(args.manifest_file, args.brc_mode, args.ignore_final_stops, args.no_fail)
    inspector.check_integrity()


if __name__ == "__main__":
    main()
