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


from typing import Any, Dict, List, Tuple
import eHive
import gzip
import json
import re
import hashlib

from functools import partial
from Bio import SeqIO
from os import path

from ensembl.brc4.runnable.seqregion_parser import SeqregionParser
from ensembl.utils.archive import open_gz_file


class SeqGroup:
    def __init__(self, sequence, identifier=None) -> None:
        self.sequence = sequence
        self.length = len(self.sequence)
        self.ids = []
        if identifier:
            self.add_id(identifier)
        self.count = len(self.ids)

    def __str__(self) -> str:
        return ", ".join(self.ids)

    def add_id(self, identifier) -> None:
        self.ids.append(identifier)
        self.count = len(self.ids)


class compare_fasta(eHive.BaseRunnable):
    def param_defaults(self):
        return {}

    def run(self) -> None:
        report = self.param_required("report")
        fasta1 = self.param_required("fasta1")
        fasta2 = self.param_required("fasta2")
        map_dna_path = self.param_required("seq_regions")
        output_dir = self.param_required("output_dir")
        species = self.param_required("species")
        name = self.param_required("comparison_name")
        accession = self.param_required("accession")

        map_dna = self.get_map(map_dna_path)
        seq1 = self.get_fasta(fasta1, map_dna)
        seq2 = self.get_fasta(fasta2, map_dna)

        (stats, diffs, seq_map) = self.compare_seqs(seq1, seq2)
        # Print mapping to a file (add report data)
        map_file = output_dir + "/" + species + "_" + name + ".map"
        self.print_map(seq_map, map_file, report, accession)

        # Print full list of results in a file
        output_file = output_dir + "/" + species + "_" + name + ".log"
        print(f"Write results in {output_file}")
        with open(output_file, "w") as out_fh:
            for line in diffs:
                out_fh.write(line + "\n")

        # Print the stats separately
        out = {"species": species, "stats": stats}
        self.dataflow(out, 2)

    def print_map(self, seq_map: dict, map_file: str, report_file: str, accession: str) -> None:
        report_parser = SeqregionParser()
        report_seq = report_parser.get_report_regions(report_file, accession)
        report = self.add_report_to_map(seq_map, report_seq)

        print(f"Write map in {map_file}")
        with open(map_file, "w") as out_fh:
            out_fh.write(json.dumps(report, sort_keys=True, indent=4))

    def add_report_to_map(self, seq_map: dict, report_seq: dict) -> List[Any]:
        accession_version = r"\.\d+$"
        report = []
        for insdc_name, old_name in seq_map.items():
            if insdc_name not in report_seq:
                raise Exception("No INSDC %s found in report" % insdc_name)
            else:
                seqr = report_seq[insdc_name]
                seqr["name"] = old_name
                seqr["EBI_seq_region_name"] = old_name
                brc4_name = insdc_name
                brc4_name = re.sub(accession_version, "", brc4_name)
                seqr["BRC4_seq_region_name"] = brc4_name
                syns = [{"source": "INSDC", "name": insdc_name}]
                seqr["synonyms"] = syns
                report.append(seqr)

        return report

    def get_map(self, map_path: str) -> dict:
        print(f"Read file {map_path}")
        data = self.get_json(map_path)

        map_dna = {}

        for seqr in data:
            name = seqr["name"]
            if "synonyms" in seqr:
                for syn in seqr["synonyms"]:
                    if syn["name"] == "INSDC":
                        map_dna[name] = syn["value"]

        return map_dna

    def get_json(self, json_path: str) -> dict:
        with open(json_path) as json_file:
            return json.load(json_file)

    def build_seq_dict(self, seqs: dict) -> dict:
        """Build a seq dict taking duplicates into account"""

        seqs_dict = dict()
        for name, seq in seqs.items():
            if seq in seqs_dict:
                seqs_dict[seq].add_id(name)
            else:
                seqs_dict[seq] = SeqGroup(seq, name)

        return seqs_dict

    def get_fasta(self, fasta_path: str, map_dna: dict) -> dict:
        print(f"Read file {fasta_path}")
        sequences = {}
        with open_gz_file(fasta_path) as fasta_fh:
            for rec in SeqIO.parse(fasta_fh, "fasta"):
                name = rec.id
                if name in map_dna:
                    name = map_dna[name]
                sequences[name] = re.sub(r"[^CGTA]", "N", str(rec.seq.upper()))
        return sequences

    def compare_seqs(self, seq1: dict, seq2: dict) -> Tuple[dict, list, dict]:
        comp = []
        accession = self.param_required("accession")
        diff = abs(len(seq1) - len(seq2))
        stats = {
            "accession": accession,
            "seq_count_1": len(seq1),
            "seq_count_2": len(seq2),
            "num_diff_seq": diff,
            "common": 0,
            "only1": 0,
            "only2": 0,
            "max_only1": 0,
            "max_only2": 0,
            "only1_200": 0,
            "only1_1000": 0,
            "only2_200": 0,
            "only2_1000": 0,
            "other_locations": 0,
            "summary": None,
            "organellar_summary": None,
            "Assembly_level_1": None,
            "Assembly_level_2": None,
        }
        value = "identical"  # variable used for summary
        org_value = "no_organelles_present"  # variable used for organellar_summary

        # Compare sequences
        seqs1 = self.build_seq_dict(seq1)
        seqs2 = self.build_seq_dict(seq2)

        # Compare number of sequences
        if len(seq1) != len(seq2):
            comp.append(f"WARNING: Different number of sequences: {len(seq1)} vs {len(seq2)}")
        else:
            comp.append(f"Same number of sequences: {len(seq1)}")

        # Sequences that are not common
        only1 = {seq: group for seq, group in seqs1.items() if not seq in seqs2}

        only2 = {seq: group for seq, group in seqs2.items() if not seq in seqs1}

        common, group_comp = self.find_common_groups(seqs1, seqs2)
        comp += group_comp

        if only1 or only2:
            value = "mismatch"

        # Gathering the organellar sequences
        report = self.param_required("report")
        report_parser = SeqregionParser()
        report_seq = report_parser.get_report_regions(report, accession)
        map_dna_path = self.param_required("seq_regions")
        seq_data = self.get_json(map_dna_path)
        org_loc = self.organellar_assembly(report_seq, seq_data)
        INSDC_assembly_level, core_assembly_level = self.assembly_level(report_seq, seq_data)

        comp.append(f"Assembly level: {INSDC_assembly_level} vs {core_assembly_level}")

        names_length = {}
        # sequences which have extra N at the end
        if only1 and only2:
            for seq_1, name1 in only1.items():
                len1 = len(seq_1)
                seq1_N = seq_1.count("N")
                for seq_2, name2 in only2.items():
                    len2 = len(seq_2)
                    seq2_N = seq_2.count("N")
                    sequence_2 = seq_2[:len1]
                    if sequence_2 == seq_1:
                        ignored_seq = seq_2[len1:]
                        N = ignored_seq.count("N")
                        if len(ignored_seq) == N:
                            comp.append(f"Please check extra Ns added in core in {name1} and {name2}")
                        else:
                            comp.append(
                                f"ALERT INSERTIONS at the end or diff assembly level {name1} and {name2}"
                            )
                    elif len1 == len2:
                        if seq2_N > seq1_N:
                            comp.append(f"Core has more Ns, check {name1} and {name2}")
                        elif seq1_N > seq2_N:
                            comp.append(f"INSDC has more Ns, check {name1} and {name2}")
                        else:
                            names_length[name1] = name2
                    else:
                        continue

        if names_length:
            length = len(names_length)
            comp.append(f"{length} sequences have the same length")
            for insdc, core in names_length.items():
                comp.append(f"INSDC: {insdc} and coredb : {core}")

        # Remove the duplicates
        for org_name in list(org_loc.keys()):
            for insdc_id, core_id in common.items():
                if org_name == core_id:
                    org_loc.pop(org_name)

        # checking for multiple entries of organellar seq
        multi_org = [name.split(".")[0] for name in org_loc.keys()]
        multi_org_acc = [j[:-1] for j in multi_org]  # similar accession
        unique_org_id = list(set(multi_org_acc))
        location = [location for location in org_loc.values()]
        unique_location = location.count("mitochondrial_chromosome")
        unique_apicoplast = location.count("apicoplast_chromosome")

        only1_id = [str(id1) for id1 in only1.values()]

        # comparing organellar sequences with common, only1 and only2
        count = 0
        for org_name, loc in org_loc.items():
            if org_name == "na":
                comp.append("MISSING accession in the report (na)")
            else:
                if org_name in common.keys():
                    count = count + 1
                    comp.append(f"{org_name} (both) in location: {loc}")
                    if count > 0:
                        org_value = "identical"
                elif org_name in only1_id:
                    count = count + 1
                    comp.append(f"{org_name} (only1) in  location: {loc}")
                    org_value = "unknown_with_organellar"
                else:
                    count = count + 1
                    comp.append(f"{org_name} (only2) in location: {loc}")
                    org_value = "unknown_with_organellar"

        # if the mistmatch is due to added organellar sequences
        if len(seqs1) > len(seqs2):
            greater_len = len(seq1)
        else:
            greater_len = len(seq2)

        diff_common = greater_len - len(common)
        diff = abs(len(only1) + len(only2))

        if diff != 0:
            if diff == count and diff_common == count:
                org_value = "organellar_present"

        if count == 0:
            org_value = "no_organelles_present"

        # checking if multiple entries of organellar sequences are present
        if len(multi_org_acc) != len(unique_org_id):
            if unique_location > 1 or unique_apicoplast > 1:
                org_value = "WARNING:Multiple_entry"

        # updating the stats
        stats["num_diff_seq"] = diff
        stats["common"] = len(common)
        stats["only1"] = len(only1)
        stats["only2"] = len(only2)
        stats["other_locations"] = count
        stats["summary"] = value
        stats["organellar_summary"] = org_value
        stats["Assembly_level_1"] = INSDC_assembly_level
        stats["Assembly_level_2"] = core_assembly_level
        print(stats)

        if only1:
            stats["max_only1"] = len(max(only1, key=lambda k: len(k)))
            # Only list sequences where the length is > 200
            mini = {seq: name for seq, name in only1.items() if len(seq) <= 200}
            maxi = {seq: name for seq, name in only1.items() if len(seq) > 200}

            if mini and len(mini) > 3000:
                comp.append(f"WARNING: Ignoring {len(mini)} sequences from 1 with length <= 200")
                only1 = maxi

        if only1:
            # Only list sequences where the length is > 1000
            mini = {seq: name for seq, name in only1.items() if len(seq) <= 1000}
            maxi = {seq: name for seq, name in only1.items() if len(seq) > 1000}
            if mini and len(mini) > 3000:
                comp.append(f"WARNING: Ignoring {len(mini)} sequences from 1 with length <= 1000")
                only1 = maxi

        if only1:
            total = sum([len(seq) for seq in only1.keys()])
            comp.append(f"WARNING: Sequences only in 1: {len(only1)} ({total})")
            only_seq1 = {name: len(seq) for seq, name in only1.items()}
            for name, length in sorted(only_seq1.items(), key=lambda x: x[1]):
                comp.append(f"\tOnly in 1: {name} ({length})")

        if only2:
            stats["max_only2"] = len(max(only2, key=lambda k: len(k)))
            # Only list sequences where the length is > 200
            mini = {seq: name for seq, name in only2.items() if len(seq) <= 200}
            maxi = {seq: name for seq, name in only2.items() if len(seq) > 200}

            if mini and len(mini) > 3000:
                comp.append(f"WARNING: Ignoring {len(mini)} sequences from 2 with length <= 200")
                only2 = maxi

        if only2:
            # Only list sequences where the length is > 1000
            mini = {seq: name for seq, name in only2.items() if len(seq) <= 1000}
            maxi = {seq: name for seq, name in only2.items() if len(seq) > 1000}

            if mini and len(mini) > 3000:
                comp.append(f"WARNING: Ignoring {len(mini)} sequences from 2 with length <= 1000")
                only2 = maxi

        if only2:
            total = sum([len(seq) for seq in only2.keys()])
            comp.append(f"WARNING: Sequences only in 2: {len(only2)} ({total})")
            only_seq2 = {name: len(seq) for seq, name in only2.items()}
            for name, length in sorted(only_seq2.items(), key=lambda x: x[1]):
                comp.append(f"\tOnly in 2: {name} ({length})")

        return (stats, comp, common)

    def find_common_groups(self, seqs1: dict, seqs2: dict) -> Tuple[dict, List[Any]]:
        print(len(seqs1))
        print(len(seqs2))
        comp = []
        common = {}
        for seq1, group1 in seqs1.items():
            if seq1 in seqs2:
                group2 = seqs2[seq1]
                # Check that the 2 groups have the same number of sequences
                if group1.count == group2.count:
                    if group1.count == 1:
                        common[group1.ids[0]] = group2.ids[0]
                    else:
                        comp.append(f"Matched 2 identical groups of sequences: {group1} and {group2}")
                        possible_id2 = " OR ".join(group2.ids)
                        for id1 in group1.ids:
                            common[id1] = possible_id2

                else:
                    comp.append(
                        f"Matched 2 different groups of sequences ({group1.count} vs {group2.count}): {group1} and {group2}"
                    )

        print(len(common))
        return common, comp

    def organellar_assembly(self, report_seq: dict, data: List[dict]) -> dict:
        org_loc = {}

        # Gathering data from the INSDC report file and storing it into a list
        for name1, details1 in report_seq.items():
            if "location" in details1:
                if details1["location"] not in (
                    "chromosome",
                    "nuclear_chromosome",
                    "linkage_group",
                ):
                    loc = details1["location"]
                    org_loc[name1] = loc

        # Gathering data from Seq_json file and storing it into a list
        for rep in data:
            for name2, details2 in rep.items():
                if "location" in name2:
                    if details2 not in (
                        "chromosome",
                        "nuclear_chromosome",
                        "linkage_group",
                    ):
                        name = rep["BRC4_seq_region_name"]
                        org_loc[name] = details2

        return org_loc

    def assembly_level(self, report_seq: dict, core_data: list) -> Tuple[str, str]:
        INSDC_assembly_level = []
        core_assembly_level = []
        core_assembly = {}
        scaffold_INSDC = 0
        chromosome_INSDC = 0
        scaffold_core = 0
        chromosome_core = 0

        for name, insdc_rep in report_seq.items():
            if insdc_rep["coord_system_level"] not in (
                "chromosome",
                "nuclear_chromosome",
            ):
                scaffold_INSDC += 1
            else:
                chromosome_INSDC += 1

        INSDC_assembly_level.extend([scaffold_INSDC, chromosome_INSDC])

        for core_details in core_data:
            name = core_details["BRC4_seq_region_name"]
            coord_system_level = core_details["coord_system_level"]
            core_assembly[name] = coord_system_level

        for name, coord_level in core_assembly.items():
            if coord_level not in ("chromosome", "nuclear_chromosome"):
                scaffold_core += 1
            else:
                chromosome_core += 1

        core_assembly_level.extend([scaffold_core, chromosome_core])

        INSDC_assembly_level = ", ".join([str(assembly) for assembly in INSDC_assembly_level])
        core_assembly_level = ", ".join([str(assembly) for assembly in core_assembly_level])

        return INSDC_assembly_level, core_assembly_level
