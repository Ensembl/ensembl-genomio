#!env python3
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


import gzip
import hashlib
import io
import json
from math import floor
from os import path
from pathlib import Path
import re
import sys

from BCBio import GFF
from Bio import SeqIO
import eHive

from ensembl.io.genomio.utils.json_utils import get_json


class integrity(eHive.BaseRunnable):
    """Check the integrity of sequence and annotation files in the genome"""

    def param_defaults(self):
        return {
            "ensembl_mode": False,
            "ignore_final_stops": False,
        }

    def run(self):
        """Load files listed in the manifest.json and check the integrity.
            Check if the files are correct by verifying the MD5 hash.
            Check if translation, functional annotation and sequence region ids
            and lengths are consistent with the information in gff.
            Compare sequence length from fasta_dna file to seq_region.json metadata.

        Args:
            manifest: Path to the manifest file.
            It contains a set of files fasta, json metadata
            and optional annotation files (gff, functional_annotation).

        Returns:
            Error if any of the above checks fail.
        """

        manifest_path = self.param_required("manifest")
        errors = []
        # load the manisfest.json
        with open(manifest_path) as manifest_file:
            manifest = json.load(manifest_file)

            # Use dir name from the manifest
            for name in manifest:
                if "file" in manifest[name]:
                    file_name = manifest[name]["file"]
                    file_name = path.join(path.dirname(manifest_path), file_name)
                    # check if the md5sum is correct
                    md5sum = manifest[name]["md5sum"]
                    errors += self.check_md5sum(file_name, md5sum)

                    manifest[name] = file_name
                else:
                    for f in manifest[name]:
                        if "file" in manifest[name][f]:
                            file_name = manifest[name][f]["file"]
                            file_name = path.join(path.dirname(manifest_path), file_name)
                            # check if the md5sum is correct
                            md5sum = manifest[name][f]["md5sum"]
                            errors += self.check_md5sum(file_name, md5sum)

                            manifest[name][f] = file_name

            # Get content from the manifest file and store it into following variables
            dna = {}
            pep = {}
            seq_regions = {}
            seq_lengths = {}
            seq_circular = {}
            gff = {}
            func_ann = {}
            agp_seqr = {}
            genome = {}

            if "gff3" in manifest:
                print("Got a gff")
                gff = self.get_gff3(manifest["gff3"])
            if "fasta_dna" in manifest:
                print("Got a fasta dna")
                # Verify if the length and id for the sequence is unique
                dna, dna_errors = self.get_fasta_lengths(manifest["fasta_dna"])
                errors += dna_errors
            if "fasta_pep" in manifest:
                print("Got a fasta pep")
                ignore_final_stops = self.param("ignore_final_stops")
                # Verify if the length and id for the sequence is unique
                pep, pep_errors = self.get_fasta_lengths(
                    manifest["fasta_pep"], ignore_final_stops=ignore_final_stops
                )
                errors += pep_errors
            if "seq_region" in manifest:
                print("Got a seq_regions")
                seq_regions = get_json(Path(manifest["seq_region"]))
                seqr_lengths = {}
                seqr_seqlevel = {}
                # Store the length as int
                for seq in seq_regions:
                    seq_lengths[seq["name"]] = int(seq["length"])
                    seq_circular[seq["name"]] = seq.get("circular", False)
                    if seq["coord_system_level"] == "contig":
                        seqr_seqlevel[seq["name"]] = int(seq["length"])
            if "functional_annotation" in manifest:
                print("Got a func_anns")
                func_ann = self.get_functional_annotation(manifest["functional_annotation"])
            if "agp" in manifest:
                print("Got agp files")
                agp_seqr = self.get_agp_seq_regions(manifest["agp"])
            if "genome" in manifest:
                print("Got a genome")
                genome = get_json(Path(manifest["genome"]))

            # Check if the accession is correct in genome.json
            if genome:
                if "assembly" in genome:
                    genome_ass = genome["assembly"]
                    if "accession" in genome_ass:
                        genome_acc = genome_ass["accession"]
                        if not re.match("GC[AF]_\d{9}(\.\d+)?", genome_acc):
                            errors += ["Genome assembly accession is wrong: '%s'" % genome_acc]

            # Check gff3
            if gff:
                # Check fasta_pep.fa integrity
                # The sequence length and id retrieved from the fasta_pep file and compared to the translated CDS id and length in the gff
                # We don't compare the peptide lengths because of seqedits
                if pep:
                    tr_errors = self.check_lengths(
                        pep, gff["translations"], "Fasta translations vs gff", special_diff=True
                    )
                    if len(tr_errors) > 0:
                        # The pseudo CDSs are included in this check
                        # Pseudo CDSs are not translated, if the pseudo translation ids are not ignored in the gff it will give an error
                        tr_errors = self.check_lengths(
                            pep,
                            gff["all_translations"],
                            "Fasta translations vs gff (include pseudo CDS)",
                            special_diff=True,
                        )
                        errors += tr_errors

                # Check functional_annotation.json integrity
                # Gene ids, translated CDS ids and translated CDSs including pseudogenes are compared to the gff
                if func_ann:
                    errors += self.check_ids(func_ann["genes"], gff["genes"], "Gene ids metadata vs gff")
                    tr_errors = self.check_ids(
                        func_ann["translations"], gff["translations"], "Translation ids metadata vs gff"
                    )
                    if len(tr_errors) > 0:
                        tr_errors = self.check_ids(
                            func_ann["translations"],
                            gff["all_translations"],
                            "Translation ids metadata vs gff (include pseudo CDS)",
                        )
                    errors += tr_errors
                    errors += self.check_ids(
                        func_ann["transposable_elements"],
                        gff["transposable_elements"],
                        "TE ids metadata vs gff",
                    )

                # Check the seq.json intregrity
                # Compare the length and id retrieved from seq.json to the gff
                if seq_regions:
                    errors += self.check_seq_region_lengths(
                        seq_lengths, gff["seq_region"], "Seq_regions metadata vs gff", seq_circular
                    )

            # Check fasta dna and seq_region integrity
            if dna and seq_regions:
                errors += self.check_seq_region_lengths(seq_lengths, dna, "seq_regions json vs dna")

            # Check agp and seq_region integrity
            if agp_seqr and seq_lengths:
                errors += self.check_seq_region_lengths(seq_lengths, agp_seqr, "seq_regions json vs agps")

        if errors:
            errors_str = "\n".join(errors)
            raise Exception("Integrity test failed for %s:\n%s" % (manifest_path, errors_str))

    def check_md5sum(self, path, md5sum):
        """Verify the integrity of the files in manifest.json.

            An MD5 hash is generated using the path provided which is then compared to the hash
            in manifest.json.

        Args:
            Path: The path for each file in the genome.
            md5sum: MD5 hash for the files.

        Returns:
            Error if the md5sum does not match.
        """

        errors = []
        with open(path, "rb") as f:
            bytes = f.read()
            readable_hash = hashlib.md5(bytes).hexdigest()
            if readable_hash != md5sum:
                errors.append("File %s has a wrong md5sum" % path)

        return errors

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
        for rec in SeqIO.parse(fasta_path, "fasta"):
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

        errors = []
        if empty_id_count > 0:
            errors.append("%d sequences with empty ids in %s" % (empty_id_count, fasta_path))
        if non_unique_count > 0:
            errors.append("%d non unique sequence ids in %s" % (non_unique_count, fasta_path))
        if contains_stop_codon > 0:
            errors.append("%d sequences with stop codons in %s" % (contains_stop_codon, fasta_path))
        return data, errors

    def get_functional_annotation(self, json_path):
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
            tes = {}

            for item in data:
                if item["object_type"] == "gene":
                    genes[item["id"]] = 1
                elif item["object_type"] == "translation":
                    translations[item["id"]] = 1
                if item["object_type"] == "transposable_element":
                    tes[item["id"]] = 1

            return {"genes": genes, "translations": translations, "transposable_elements": tes}

    def get_gff3(self, gff3_path):
        # Load the gff file
        if gff3_path.endswith(".gz"):
            with io.TextIOWrapper(gzip.open(gff3_path, "r")) as gff3_handle:
                return self.parse_gff3(gff3_handle)
        else:
            with open(gff3_path, "r") as gff3_handle:
                return self.parse_gff3(gff3_handle)

    def parse_gff3(self, gff3_handle):
        """A GFF parser is used to retrieve information in the GFF file such as
           gene and CDS ids and their corresponding lengths.

        Args:
            gff3_handle: Path to gff3 file.

        Returns:
            dict containing sequence ids, gene ids, transcript ids and translation ids
            are stored with their corresponding lengths.
        """

        ensembl_mode = self.param("ensembl_mode")
        seqs = {}
        genes = {}
        peps = {}
        all_peps = {}
        tes = {}

        gff = GFF.parse(gff3_handle)
        for seq in gff:
            seqs[seq.id] = len(seq.seq)

            for feat in seq.features:
                feat_length = abs(feat.location.end - feat.location.start)
                # Store gene id and length
                if feat.type in ["gene", "ncRNA_gene", "pseudogene"]:
                    gene_id = feat.id
                    if ensembl_mode:
                        gene_id = feat_length
                    genes[gene_id] = abs(feat.location.end - feat.location.start)
                    # Get CDS id and length
                    for feat2 in feat.sub_features:
                        if feat2.type in ("mRNA", "pseudogenic_transcript"):
                            length = {}
                            for feat3 in feat2.sub_features:
                                if feat3.type == "CDS":
                                    pep_id = feat3.id
                                    if ensembl_mode:
                                        pep_id = pep_id.replace("CDS:", "")
                                    if pep_id not in length:
                                        length[pep_id] = 0
                                    length[pep_id] += abs(feat3.location.end - feat3.location.start)
                            for pep_id in length:
                                # Store length for translations, add pseudo translations separately
                                pep_length = floor(length[pep_id] / 3) - 1
                                if feat.type != "pseudogene":
                                    peps[pep_id] = pep_length
                                all_peps[pep_id] = pep_length
                if feat.type == "transposable_element":
                    tes[feat.id] = feat_length

        stats = {
            "seq_region": seqs,
            "genes": genes,
            "translations": peps,
            "all_translations": all_peps,
            "transposable_elements": tes,
        }
        return stats

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
                        asm_start,
                        asm_end,
                        asm_part,
                        typ,
                        cmp_id,
                        cmp_start,
                        cmp_end,
                        cmp_strand,
                    ) = line.split("\t")
                    # Ignore WGS contig
                    if typ != "W":
                        continue

                    # Assembled seq length
                    if not asm_id in seqr or seqr[asm_id] < int(asm_end):
                        seqr[asm_id] = int(asm_end)

                    # Composite seq length
                    if not cmp_id in seqr or seqr[cmp_id] < int(cmp_end):
                        seqr[cmp_id] = int(cmp_end)

        return seqr

    def check_ids(self, list1, list2, name):
        """Compare the ids in list1 and list2.

        Args:
            list1: dict containing sequence ids retrieved from functional.json.
            list2: dict containing length and id in the retrieved from the gff.
            name:  string

        Return:
            Error if the ids in functional.json and gff do not match.
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
            print("%d common elements in %s" % (len(common), name))
        if only1:
            errors.append("%d only in first list in %s (first: %s)" % (len(only1), name, only1[0]))
        if only2:
            errors.append("%d only in second list in %s (first: %s)" % (len(only2), name, only2[0]))

        return errors

    def check_lengths(self, list1, list2, name, allowed_len_diff=None, special_diff=False):
        """Check the difference in ids and length between list1 and list2.
            There are a few special cases here where we allow a certain asymmetry
            by changing the values of the arguments.

        Args:
            list1: dict containing length and id of the sequence from fasta files.
            list2: dict containing length and id in the retrieved from the gff.
            name:  string

        allowed_len_diff : set as None when we do not want to accept any difference in length between list1 and list2.
            The value here can be changed based on how much difference in sequence length we are wanting to accept.

        special_diff: set as False when no special length difference is expected between the lists.
                    This can be changed if we want to report common sequences with 1 BP difference.

        Returns:
            Error if there is a difference in length or ids between the lists.
        """

        # check list diffferences, checks if abs(values diff) < allowed_len_diff

        set1 = frozenset(list1)
        set2 = frozenset(list2)
        list1_2 = list(set1 - set2)
        list2_1 = list(set2 - set1)

        errors = []
        if len(list1_2) > 0:
            errors.append("%s: %d from the first list only (i.e. %s)" % (name, len(list1_2), list1_2[0]))
        if len(list2_1) > 0:
            errors.append("%s: %d from the second list only (i.e. %s)" % (name, len(list2_1), list2_1[0]))

        common_len = 0
        if allowed_len_diff is None:
            common_len = len(set1 & set2)
        else:
            # check for the sequence length difference
            diff_len_list = []
            diff_len_special_list = []
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
                    _dlist.append("%s: %d vs %d" % (e, list1[e], list2[e]))
            if diff_len_special_list:
                errors.append(
                    "%d common elements with one BP/AA length diff for %s (e.g. %s)"
                    % (len(diff_len_special_list), name, diff_len_special_list[0])
                )
            if diff_len_list:
                errors.append(
                    "%d common elements with length diff for %s (e.g. %s)"
                    % (len(diff_len_list), name, diff_len_list[0])
                )
        if common_len > 0:
            print("%d common elements between lists for %s" % (common_len, name), file=sys.stderr)

        return errors

    def check_seq_region_lengths(self, seqrs, feats, name, circular=None):
        """Check the integrity of seq_region.json file by comparing the length of the sequence
            to fasta files and the gff.

            Seq_region file is in json format containing the metadata of the sequence.
            It contains sequence id, length, location and the synonyms for the sequence name from different sources.

        Args:
            seqs: Sequence name and length retrieved from seq_region.json file.
            feats: Sequence name and length retrieved from the fasta and gff file.
            name: String

        Returns:
            Error if there are common sequences with difference in ids
            and if the sequences are not consistent in the files.
        """

        only_seqr = []
        only_feat = []

        common = []
        diff = []
        diff_list = []

        # not an error on circular for gff features
        diff_circular = []
        diff_circular_list = []

        for seq_id in seqrs:
            if seq_id in feats:
                # Check that feature is within the seq_region length
                if feats[seq_id] > seqrs[seq_id]:
                    if circular is None or not circular.get(seq_id, False):
                        diff.append(seq_id)
                        diff_list.append("%s: %d vs %d" % (seq_id, seqrs[seq_id], feats[seq_id]))
                    else:
                        diff_circular.append(seq_id)
                        diff_circular_list.append("%s: %d vs %d" % (seq_id, seqrs[seq_id], feats[seq_id]))
                else:
                    common.append(seq_id)
            else:
                only_seqr.append(seq_id)

        for seq_id in feats:
            if seq_id not in common and seq_id not in diff and seq_id not in diff_circular:
                only_feat.append(seq_id)

        if common:
            print("%d common elements in %s" % (len(common), name))
        if diff_circular:
            print(
                "%d differences for circular elements in %s (e.g. %s)"
                % (len(diff_circular), name, diff_circular_list[0])
            )

        errors = []
        if diff:
            errors.append(
                "%d common elements with higher length in %s (e.g. %s)" % (len(diff), name, diff_list[0])
            )
        if only_seqr:
            # Not an error!
            print("%d only in seq_region list in %s (first: %s)" % (len(only_seqr), name, only_seqr[0]))
        if only_feat:
            errors.append("%d only in second list in %s (first: %s)" % (len(only_feat), name, only_feat[0]))

        return errors
