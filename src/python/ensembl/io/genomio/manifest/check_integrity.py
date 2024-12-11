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
"""Compare the genomic data in a DNA FASTA file, seq_region JSON, gene models GFF3 and peptide FASTA
to ensure their contents are in sync.
"""

__all__ = ["IntegrityTool"]

import logging
from pathlib import Path
import re
from typing import Any

import ensembl.io.genomio
from ensembl.io.genomio.manifest.manifest_stats import InvalidIntegrityError, ManifestStats
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


class IntegrityTool:
    """Check the integrity of sequence and annotation files in the genome"""

    def __init__(
        self,
        manifest_file: Path,
        ignore_final_stops: bool = False,
        no_fail: bool = False,
    ) -> None:
        self.manifest = ManifestStats(manifest_file)
        self.ignore_final_stops = False
        self.set_ignore_final_stops(ignore_final_stops)
        self.errors: list[str] = []
        self.no_fail = no_fail

    def add_errors(self, errors: list[str] | str) -> None:
        """Store the given errors (list or single string) in the list of all errors."""
        if isinstance(errors, str):
            self.errors.append(errors)
        else:
            self.errors += errors

    def check_integrity(self) -> None:
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

            self.check_seq_region_lengths(
                seq_lengths, gff_seq_regions, "seq_regions JSON vs GFF3 lengths", seq_circular
            )

        self.check_seq_region_lengths(seq_lengths, dna, "seq_regions JSON vs DNA lengths")
        self.check_seq_region_lengths(seq_lengths, agp_seqr, "seq_regions JSON vs AGPs lengths")

        if self.errors:
            errors_str = "\n".join(self.errors)
            message = f"Integrity test failed:\n{errors_str}"
            if self.no_fail:
                print(message)
            else:
                raise InvalidIntegrityError(message)

    def set_ignore_final_stops(self, ignore_final_stops: bool) -> None:
        """Set ignore_final_stops (when calculating peptide length) for this tool and the manifest."""
        self.ignore_final_stops = ignore_final_stops
        self.manifest.ignore_final_stops = ignore_final_stops

    def _check_genome(self, genome: dict[str, Any]) -> None:
        """Check if the accession is correct in genome.json."""
        genome_accession = genome.get("assembly", {}).get("accession", "")
        if not genome_accession:
            return
        if not re.match(r"GC[AF]_\d{9}(\.\d+)?", genome_accession):
            self.add_errors(f"Genome assembly accession is wrong: '{genome_accession}'")

    def check_ids(self, list1: dict[str, Any], list2: dict[str, Any], name: str) -> list[str]:
        """Compare the ids in list1 and list2.

        Args:
            list1: Sequence IDs retrieved from `functional.json`.
            list2: Sequence IDs retrieved from the GFF3 file.
            name: Source name.

        Return:
            List of message errors of sequence IDs found only in one of the lists provided.
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

    def check_lengths(
        self,
        list1: dict[str, int],
        list2: dict[str, int],
        name: str,
        *,
        allowed_len_diff: int | None = None,
        special_diff: bool = False,
    ) -> list[str]:
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
            diff_len_list: list[str] = []
            diff_len_special_list: list[str] = []
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
        seqrs: dict[str, Any] | None,
        feats: dict[str, Any] | None,
        name: str,
        circular: dict[str, Any] | None = None,
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
        if not seqrs or not feats:
            return
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
        self, seqrs: dict[str, Any], feats: dict[str, Any], circular: dict[str, Any] | None = None
    ) -> dict[str, list[str]]:
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
        comp: dict[str, list[str]] = {
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
    parser.add_argument(
        "--ignore_final_stops", action="store_true", help="Ignore final stop when calculating peptide length"
    )
    parser.add_argument(
        "--no_fail", action="store_true", help="In case of errors, don't fail but print errors to stdout."
    )
    parser.add_argument("--version", action="version", version=ensembl.io.genomio.__version__)
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    inspector = IntegrityTool(args.manifest_file, args.ignore_final_stops, args.no_fail)
    inspector.check_integrity()


if __name__ == "__main__":
    main()
