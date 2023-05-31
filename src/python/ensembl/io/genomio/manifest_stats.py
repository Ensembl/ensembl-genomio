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
"""Compute stats from the current genome files associated with the manifest.
"""

import gzip
import json
from pathlib import Path
from shutil import which
from statistics import mean
import subprocess
from typing import Dict, List, Optional, Set, TextIO, Union

import argschema
from BCBio import GFF


class BiotypeCounter:
    """A counter for a given biotype, given a list of features."""

    def __init__(self, count: int = 0, ids: Optional[Set[str]] = None, example: Optional[str] = None) -> None:
        self.count: int = count
        if ids is None:
            ids = set()
        self.ids: Set[str] = ids
        if example is None:
            example = ""
        self.example: str = example

    def add_id(self, feature_id: str) -> None:
        """Add a feature to the counter.

        Args:
            feature_id (str): Feature id to add.
        """
        self.count += 1
        self.ids.add(feature_id)

    def unique_count(self) -> int:
        """Total number feature ids added to the counter so far.

        Returns:
            int: number of features in the counter.
        """
        return len(self.ids)


class StatsError(Exception):
    """Raised when stats could not be computed."""


class manifest_stats:
    """Representation of the statistics of the set of files listed in the manifest file provided."""

    def __init__(self, manifest_dir: str, accession: Optional[str], datasets_bin: Optional[str]):
        self.manifest = f"{manifest_dir}/manifest.json"
        self.accession: Optional[str] = accession
        self.error = False
        if datasets_bin is None:
            datasets_bin = "datasets"
        self.datasets_bin = datasets_bin
        self.manifest_parent = manifest_dir
        self.check_ncbi = False

    def run(self):
        """Compute stats in the files and output a stats.txt file in the same folder.

        Raises:
            StatsError: Could not compute some stats.
        """
        manifest = self.get_manifest()

        stats = []
        if self.accession is not None:
            stats.append(self.accession)

        if not self.error:
            if "gff3" in manifest:
                stats += self.get_gff3_stats(Path(manifest["gff3"]))
            if "seq_region" in manifest:
                stats += self.get_seq_region_stats(Path(manifest["seq_region"]))

        stats_path = f"{self.manifest_parent}/stats.txt"
        print(stats_path)
        with open(stats_path, "w") as stats_out:
            stats_out.write("\n".join(stats))

        # Flow out if errors in stats comparison
        if self.error:
            raise StatsError(f"Stats count errors, check the file {stats_path}")

    def get_manifest(self) -> Dict:
        """Get the files metadata from the manifest json file.

        Returns:
            Dict: A representation of the manifest json data.
        """
        with open(self.manifest) as f_json:
            manifest = json.load(f_json)
            manifest_root = self.manifest_parent

        # Use dir name from the manifest
        for name in manifest:
            if "file" in manifest[name]:
                file_name = manifest[name]["file"]
                file_name = f"{manifest_root}/{file_name}"
                manifest[name] = file_name
            else:
                for f in manifest[name]:
                    if "file" in manifest[name][f]:
                        file_name = manifest[name][f]["file"]
                        file_name = manifest_root, file_name
                        manifest[name][f] = file_name

        return manifest

    def get_seq_region_stats(self, seq_region_path: Path) -> List[str]:
        """Compute stats from the seq_region json file.

        Args:
            seq_region_path (Path): the seq_region json file.

        Returns:
            List[str]: Stats from the seq_regions.
        """
        json_file = open(seq_region_path, "r")
        seq_regions = json.load(json_file)

        # Get basic data
        coord_systems: Dict[str, List[int]] = {}
        circular = 0
        locations = []
        codon_tables = []
        for seqr in seq_regions:
            # Get readable seq_region name
            genbank = "synonyms" in seqr and [x for x in seqr["synonyms"] if x["source"] == "GenBank"]
            seqr_name = genbank and genbank[0]["name"] or seqr["name"]

            coord_level = seqr["coord_system_level"]
            if coord_level not in coord_systems:
                coord_systems[coord_level] = []
            coord_systems[coord_level].append(seqr["length"])

            if "circular" in seqr:
                circular += 1
            if "codon_table" in seqr:
                codon_tables.append(f"{seqr_name} = {seqr['codon_table']}")
            if "location" in seqr:
                locations.append(f"{seqr_name} = {seqr['location']}")

        # Stats
        stats: List[str] = []
        stats.append(seq_region_path.name)
        stats.append(f"Total coord_systems {len(coord_systems)}")
        for coord_name, lengths in coord_systems.items():
            stats.append(f"\nCoord_system: {coord_name}")

            stat_counts: Dict[str, Union[int, float]] = dict()
            stat_counts["Number of sequences"] = len(lengths)
            stat_counts["Sequence length sum"] = sum(lengths)
            stat_counts["Sequence length minimum"] = min(lengths)
            stat_counts["Sequence length mean"] = mean(lengths)
            stat_counts["Sequence length maximum"] = max(lengths)

            for name, count in stat_counts.items():
                stats.append(f"{count: 9d}\t{name}")

        # Special
        if circular or locations:
            stats.append("\nSpecial")
            if circular:
                stats.append(f"{circular: 9d}\tcircular sequences")
            if locations:
                stats.append(f"{len(locations): 9d} sequences with location")
                for loc in locations:
                    stats.append(f"\t\t\t%s{loc}")
            if codon_tables:
                stats.append(f"{len(codon_tables): 9d} sequences with codon_table")
                for table in codon_tables:
                    stats.append(f"\t\t\t{table}")

        stats.append("\n")

        return stats

    def get_gff3_stats(self, gff3_path: Path) -> List[str]:
        """Compute stats from the gene models in the GFF3 file.

        Args:
            gff3_path (Path): the GFF3 file.

        Returns:
            List[str]: Stats from the gene models.
        """
        stats: List[str] = []
        stats.append(gff3_path.name)
        if gff3_path.name.endswith(".gz"):
            with gzip.open(gff3_path, "rt") as gff3_handle:
                stats += self.parse_gff3(gff3_handle)
        else:
            with gff3_path.open("r") as gff3_handle:
                stats += self.parse_gff3(gff3_handle)
        stats.append("\n")

        return stats

    def parse_gff3(self, gff3_handle: TextIO) -> List:
        """Extract the gene models from the GFF3 file and compute stats.

        Args:
            gff3_handle (TextIO): the GFF3 file.

        Returns:
            List: Stats from the gene model.
        """
        biotypes: Dict[str, BiotypeCounter] = {}

        for rec in GFF.parse(gff3_handle):
            for feat1 in rec.features:
                # Check if the gene contains proteins (CDSs),
                # and keep a count of all hierarchies (e.g. gene-mRNA-CDS)
                is_protein = False
                for feat2 in feat1.sub_features:
                    if feat2.type == "mRNA":
                        types2 = {f.type for f in feat2.sub_features}
                        if "CDS" in types2:
                            is_protein = True
                    manifest_stats.increment_biotype(biotypes, feat2.id, f"{feat1.type}-{feat2.type}")
                    for feat3 in feat2.sub_features:
                        if feat3.type == "exon":
                            continue
                        manifest_stats.increment_biotype(
                            biotypes, feat3.id, f"{feat1.type}-{feat2.type}-{feat3.type}"
                        )

                # Main categories counts
                if feat1.type == "pseudogene":
                    manifest_stats.increment_biotype(biotypes, feat1.id, "pseudogene")
                elif is_protein:
                    manifest_stats.increment_biotype(biotypes, feat1.id, f"PROT_{feat1.type}")
                else:
                    # Special case, undefined gene-transcript
                    if (
                        feat1.type == "gene"
                        and feat1.sub_features
                        and feat1.sub_features[0].type == "transcript"
                    ):
                        manifest_stats.increment_biotype(biotypes, feat1.id, "OTHER")
                    else:
                        manifest_stats.increment_biotype(biotypes, feat1.id, f"NONPROT_{feat1.type}")

                # Total
                if feat1.type in ("gene", "pseudogene"):
                    manifest_stats.increment_biotype(biotypes, feat1.id, "ALL_GENES")

        # Order
        sorted_biotypes = dict()
        for name in sorted(biotypes.keys()):
            data: BiotypeCounter = biotypes[name]
            sorted_biotypes[name] = data

        stats = [
            f"{data.unique_count():>9}\t{biotype:<20}\tID = {data.example}"
            for (biotype, data) in sorted_biotypes.items()
        ]

        # Check against NCBI stats
        stats += self.check_ncbi_stats(biotypes)

        return stats

    def check_ncbi_stats(self, biotypes: Dict[str, BiotypeCounter]) -> List[str]:
        """Use the dataset tool from NCBI to get stats and compare with what we have"""
        stats: List[str] = []
        if not self.check_ncbi:
            return stats

        if self.accession is None:
            return stats

        accession: str = self.accession

        datasets_bin = self.datasets_bin
        if not which(datasets_bin):
            return stats

        # Get the dataset summary from NCBI
        command = [datasets_bin, "summary", "genome", "accession", accession]
        result_out = subprocess.run(command, stdout=subprocess.PIPE, check=True)
        result = json.loads(result_out.stdout)

        # Get stats
        if "reports" in result:
            genome = result["reports"][0]
            if "annotation_info" in genome and "stats" in genome["annotation_info"]:
                ncbi_stats = genome["annotation_info"]["stats"]

                if "gene_counts" in ncbi_stats:
                    counts = ncbi_stats["gene_counts"]
                    stats = self.compare_ncbi_counts(biotypes, counts)
        return stats

    def compare_ncbi_counts(self, biotypes: Dict[str, BiotypeCounter], ncbi: Dict) -> List[str]:
        """Compare specific gene stats from NCBI"""
        stats: List[str] = []

        maps = [
            ["total", "ALL_GENES"],
            ["protein_coding", "PROT_gene"],
            ["pseudogene", "pseudogene"],
            ["non_coding", "NONPROT_gene"],
            ["other", "OTHER"],
        ]

        for count_map in maps:
            ncbi_name, prep_name = count_map
            ncbi_count = ncbi.get(ncbi_name, 0)
            preped: Optional[BiotypeCounter] = biotypes.get(prep_name)
            prep_count = 0
            if preped is not None:
                prep_count = preped.count

            if prep_count != ncbi_count:
                diff = prep_count - ncbi_count
                stats.append(f"DIFF gene count for {count_map}: {prep_count} - {ncbi_count} = {diff}")
                self.error = True
            else:
                stats.append(f"Same count for {count_map}: {prep_count}")

        return stats

    @staticmethod
    def increment_biotype(biotypes: Dict[str, BiotypeCounter], feature_id: str, feature_biotype: str) -> None:
        """Add the feature to their respective biotype counter.

        Args:
            biotypes (Dict[str, BiotypeCounter]): All current biotypes, with their counter.
            feature_id (str): Feature id to be counted.
            feature_biotype (str): The biotype of the feature.
        """
        if feature_biotype not in biotypes:
            biotypes[feature_biotype] = BiotypeCounter(example=feature_id)
        biotypes[feature_biotype].add_id(feature_id)


class InputSchema(argschema.ArgSchema):  # need more metadata/'True' requirements
    """Input arguments expected by this script."""

    manifest_dir = argschema.fields.String(required=True, metadata={"description": "Manifest file path"})
    accession = argschema.fields.String(
        metadata={"description": "Sequence accession ID to compare stats with NCBI"}
    )
    datasets_bin = argschema.fields.String(metadata={"description": "Datasets bin status"})


def main():
    """Main entrypoint."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    mstats = manifest_stats(mod.args["manifest_dir"], mod.args.get("accession"), mod.args.get("datasets_bin"))
    if mod.args.get("accession"):
        mstats.check_ncbi = True
    mstats.run()


if __name__ == "__main__":
    main()
