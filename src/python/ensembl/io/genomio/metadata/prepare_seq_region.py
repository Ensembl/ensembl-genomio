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
"""TODO"""

__all__ = [
    "SeqRegion",
    "LOCATION_CODON",
    "MOLECULE_LOCATION",
    "SYNONYM_MAP",
    "SYNONYM_RESOURCES",
    "add_insdc_seq_region_name",
    "add_mitochondrial_codon_table",
    "add_translation_table",
    "exclude_seq_regions",
    "get_codon_table",
    "get_gbff_seq_regions",
    "get_genbank_id",
    "get_organelle",
    "get_report_regions",
    "make_seq_region",
    "merge_seq_regions",
    "prepare_seq_region_metadata",
    "report_to_csv",
]

import ast
import csv
import gzip
from os import PathLike
from pathlib import Path
import re
import requests
from typing import Any, Dict, List, Optional, Tuple

import argschema
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from marshmallow import post_load

from ensembl.io.genomio.utils import get_json, print_json


# Definition of SeqRegion type
SeqRegion = Dict[str, Any]


SYNONYM_RESOURCES = ["GenBank", "RefSeq", "INSDC"]
SYNONYM_MAP = {
    "Assigned-Molecule": "INSDC",
    "GenBank-Accn": "GenBank",
    "RefSeq-Accn": "RefSeq",
    "Sequence-Name": "INSDC_submitted_name",
}
MOLECULE_LOCATION = {
    "apicoplast": "apicoplast_chromosome",
    "chromosome": "nuclear_chromosome",
    "kinetoplast": "kinetoplast_chromosome",
    "linkage group": "linkage_group",
    "mitochondrion": "mitochondrial_chromosome",
    "plasmid": "plasmid",
}
LOCATION_CODON = {"apicoplast_chromosome": 4}


def exclude_seq_regions(seq_regions: List[SeqRegion], to_exclude: List[str]) -> List[SeqRegion]:
    """Returns the list of sequence regions with the ones from the exclusion list removed.

    Args:
        seq_regions: Sequence regions.
        to_exclude: Sequence region names to exclude.

    """
    filtered_seq_regions = []
    for seqr in seq_regions:
        if ("name" in seqr) and (seqr["name"] in to_exclude):
            print(f'Remove seq_region {seqr["name"]}')
        else:
            filtered_seq_regions.append(seqr)
    return filtered_seq_regions


def add_translation_table(
    seq_regions: List[SeqRegion], location_codon: Dict[str, int] = LOCATION_CODON
) -> None:
    """Adds the translation codon table to each sequence region (when missing) based on its location.

    Args:
        seq_regions: Sequence regions.
        location_codon: Map of known codon tables for known locations.

    """
    for seqr in seq_regions:
        # Do not overwrite any existing codon table
        if ("codon_table" not in seqr) and ("location" in seqr) and (seqr["location"] in location_codon):
            seqr["codon_table"] = location_codon[seqr["location"]]


def add_insdc_seq_region_name(
    seq_regions: List[SeqRegion], synonym_sources: List[str] = SYNONYM_RESOURCES
) -> List[SeqRegion]:
    """Returns the list of sequence regions with their corresponding INSDC sequence region names.

    "BRC4_seq_region_name" and "EBI_seq_region_name" fields are added to each sequence region: the
    former will contain the corresponding INSDC name whilst the latter will contain the current name.

    Args:
        seq_regions: Sequence regions.
        synonym_sources: Synonym sources to use for the BRC4 name, in order of preference.

    Raises:
        Exception: If no synonym name is found for a sequence region.

    """
    new_seq_regions = []
    for seqr in seq_regions:
        names = {synonym["source"]: synonym["name"] for synonym in seqr.get("synonyms", [])}
        # Choose the synonym to use as the BRC name, using the first valid name from the source list
        for source_name in synonym_sources:
            if source_name in names:
                brc_name = names[source_name]
                break
        else:
            raise Exception(f'Cannot set BRC4 sequence region name for {seqr["name"]}')
        brc_name = brc_name.partition(".")[0]
        seqr["BRC4_seq_region_name"] = brc_name
        seqr["EBI_seq_region_name"] = seqr["name"]
        new_seq_regions.append(seqr)
    return new_seq_regions


def add_mitochondrial_codon_table(seq_regions: List[SeqRegion], taxon_id: int) -> None:
    """Adds the mitochondrial codon table to each sequence region (when missing) based on the taxon ID.

    If no mitochondrial genetic code can be found for the given taxon ID nothing will be changed.

    Args:
        seq_regions: Sequence regions.
        taxon_id: The species taxon ID.

    """
    url = f"https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{str(taxon_id)}"
    response = requests.get(url, headers={"Content-Type": "application/json"})
    decoded = response.json()
    if "mitochondrialGeneticCode" not in decoded:
        print(f"No mitochondria genetic code found for taxon {taxon_id}")
    else:
        genetic_code = int(decoded["mitochondrialGeneticCode"])
        for seqr in seq_regions:
            location = seqr.get("location", "")
            # Do not overwrite any existing mitochondrial codon table
            if (location == "mitochondrial_chromosome") and ("codon_table" not in seqr):
                seqr["codon_table"] = genetic_code


def merge_seq_regions(left: Dict[str, SeqRegion] = {}, right: Dict[str, SeqRegion] = {}) -> List[SeqRegion]:
    """Merges sequence regions from different sources.

    When combining two regions matching the same key, the "right" seq region data will take precedence
    over the "left".

    Args:
        left: Dictionary of sequence regions with names as keys.
        right: Dictionary of sequence regions with names as keys.

    Returns:
        A list of merged sequence regions, sorted by their name.

    """
    # Get all the names
    all_names = set(left).union(right)
    # Create the seq_regions, merge if needed
    seq_regions = []
    for name in all_names:
        left_seqr = left.get(name, {})
        right_seqr = right.get(name, {})
        if left_seqr and right_seqr:
            final_seqr = {**left_seqr, **right_seqr}
        elif left_seqr:
            final_seqr = left_seqr
        else:
            final_seqr = right_seqr
        seq_regions.append(final_seqr)
    seq_regions.sort(key=lambda x: x["name"])
    return seq_regions


def get_gbff_seq_regions(gbff_path: PathLike) -> Dict[str, SeqRegion]:
    """Returns the sequence regions found in the GBFF file (if any).

    Args:
        gbff_path: Path to GBFF file.

    Returns:
        A dict of SeqRegions, with their name as the key.

    """
    seq_regions = {}
    _open = str(gbff_path).endswith(".gz") and gzip.open or open
    with _open(gbff_path, "rt") as gbff_file:
        for record in SeqIO.parse(gbff_file, "genbank"):
            seqr: SeqRegion = {"length": len(record.seq)}
            # Is the seq_region circular?
            annotations = record.annotations
            if ("topology" in annotations) and (annotations["topology"] == "circular"):
                seqr["circular"] = True
            # Is there a genetic code defined?
            codon_table = get_codon_table(record)
            if codon_table != None:
                seqr["codon_table"] = codon_table
            # Is it an organelle?
            location = get_organelle(record)
            if location != None:
                seqr["location"] = location
            # Is there a comment stating the Genbank record this is based on?
            genbank_id = get_genbank_id(record)
            if genbank_id != None:
                seqr["synonyms"] = [{"source": "INSDC", "name": genbank_id}]
            # Store the seq_region
            seq_regions[record.id] = seqr
    return seq_regions


def get_genbank_id(record: SeqRecord) -> Optional[str]:
    """Returns the GenBank accession from a given sequence record (if present).

    Only useful for RefSeq sequence records, where the GenBank accession is stored in a comment.

    Args:
        record: Sequence record.

    """
    genbank_id = None
    if "comment" in record.annotations:
        comment = record.annotations["comment"]
        comment = re.sub(r"[ \n\r]+", " ", comment)
        match = re.search(r"The reference sequence was derived from ([^\.]+)\.", comment)
        if match:
            genbank_id = match.group(1)
    return genbank_id


def get_codon_table(record: SeqRecord) -> Optional[int]:
    """Returns the codon table number from a given a GenBank sequence record (if present).

    Args:
        record: GenBank sequence record.

    """
    table_number = None
    for feat in record.features:
        if "transl_table" in feat.qualifiers:
            table_number = int(feat.qualifiers["transl_table"][0])
            break
    return table_number


def get_organelle(record: SeqRecord, molecule_location: Dict = MOLECULE_LOCATION) -> Optional[str]:
    """Returns the organelle location from the given GenBank record (if present).

    Args:
        record: GenBank sequence record.
        molecule_location: Map of sequence type to SO location.

    Raises:
        KeyError: If the location is not part of the controlled vocabulary.

    """
    location = None
    for feat in record.features:
        if "organelle" in feat.qualifiers:
            organelle = str(feat.qualifiers["organelle"][0])
            if not organelle:
                break
            # Remove plastid prefix
            with_prefix = re.match(r"^(plastid|mitochondrion):(.+)$", organelle)
            if with_prefix:
                organelle = with_prefix[2]
            # Get controlled name
            try:
                location = molecule_location[organelle]
            except KeyError:
                raise KeyError(f"Unrecognized sequence location: {organelle}")
            break
    return location


def get_report_regions(report_path: PathLike, is_refseq: bool) -> Dict[str, SeqRegion]:
    """Returns the sequence regions found in the report file.

    Args:
        report_path: Path to the sequence regions report file.
        is_refseq: True if the source of the report is RefSeq, false if INSDC.

    Returns:
        A dict of SeqRegions, with their name as the key.

    """
    # Get the report in a CSV format
    report_csv = report_to_csv(report_path)[0]
    # Feed the csv string to the CSV reader
    reader = csv.DictReader(report_csv.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE)
    # Create the seq_regions
    seq_regions = {}
    for row in reader:
        seq_region = make_seq_region(row, is_refseq)
        if seq_region:
            name = seq_region["name"]
            seq_regions[name] = seq_region
    return seq_regions


def make_seq_region(
    data: Dict,
    is_refseq: bool,
    synonym_map: Dict[str, str] = SYNONYM_MAP,
    molecule_location: Dict[str, str] = MOLECULE_LOCATION,
) -> SeqRegion:
    """Returns a sequence region from the information provided.

    An empty sequence region will be returned if not accession information is found.

    Args:
        data: a dict from the report representing one line, where the key is the column name.
        is_refseq: True if the source is RefSeq, false if INSDC.
        synonym_map: Map of INSDC report column names to sequence region field names.
        molecule_location: Map of sequence type to SO location.

    Raises:
        KeyError: If the sequence location is not recognised.
        Exception: if the sequence role is not recognised.

    """
    seq_region = {}
    # Set accession as the sequence region name
    src = "RefSeq" if is_refseq else "GenBank"
    accession_id = data.get(f"{src}-Accn", "")
    if accession_id and (accession_id != "na"):
        seq_region["name"] = accession_id
    else:
        print(f'No {src} accession ID found for {data["Sequence-Name"]}')
        return {}
    # Add synonyms
    synonyms = []
    for field, source in synonym_map.items():
        if (field in data) and (data[field].casefold() != "na"):
            synonym = {"source": source, "name": data[field]}
            synonyms.append(synonym)
    if synonyms:
        synonyms.sort(key=lambda x: x["source"])
        seq_region["synonyms"] = synonyms
    # Add sequence length
    field = "Sequence-Length"
    if (field in data) and (data[field].casefold() != "na"):
        seq_region["length"] = int(data[field])
    # Add coordinate system and location
    seq_role = data["Sequence-Role"]
    # Scaffold?
    if seq_role in ("unplaced-scaffold", "unlocalized-scaffold"):
        seq_region["coord_system_level"] = "scaffold"
    # Chromosome? Check location
    elif seq_role == "assembled-molecule":
        seq_region["coord_system_level"] = "chromosome"
        location = data["Assigned-Molecule-Location/Type"].lower()
        # Get location metadata
        try:
            seq_region["location"] = molecule_location[location]
        except KeyError:
            raise KeyError(f"Unrecognized sequence location: {location}")
    else:
        raise Exception(f"Unrecognized sequence role: {seq_role}")
    return seq_region


def report_to_csv(report_path: PathLike) -> Tuple[str, Dict]:
    """Returns an assembly report as a CSV string.

    Args:
        report_path: path to a seq_region file from INSDC/RefSeq

    Returns:
        The data as a string in CSV format, and the head metadata as a dictionary.

    """
    _open = str(report_path).endswith(".gz") and gzip.open or open
    with _open(report_path, "rt") as report:
        data = ""
        metadata = {}
        last_head = ""
        for line in report:
            # Ignore header
            if line.startswith("#"):
                # Get metadata values if possible
                match = re.search("# (.+?): (.+?)$", line)
                if match:
                    metadata[match.group(1)] = match.group(2)
                last_head = line
            else:
                if last_head:
                    data += last_head[2:].strip() + "\n"
                    last_head = ""
                data += line
        return data, metadata


def prepare_seq_region_metadata(
    genome_file: PathLike,
    report_file: PathLike,
    gbff_file: Optional[PathLike],
    dst_dir: PathLike,
    brc4_mode: bool = True,
    to_exclude: Optional[List[str]] = None,
) -> None:
    """Prepares the sequence region metadata found in the INSDC/RefSeq report and GBFF files.

    The sequence region information is loaded from both sources and combined. Elements are added/excluded
    as requested, and the final sequence region metadata is dumped in a JSON file that follows the schema
    defined in "schemas/seq_region_schema.json".

    Args:
        genome_file: Genome metadata JSON file path.
        report_file: INSDC/RefSeq sequences report file path to parse.
        gbff_file: INSDC/RefSeq GBFF file path to parse.
        dst_dir: Output folder for the processed sequence regions JSON file.
        brc4_mode: Include INSDC sequence region names?
        to_exclude: Sequence region names to exclude.

    """
    genome_data = get_json(genome_file)
    dst_dir = Path(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)
    # Final file name
    metadata_type = "seq_region"
    metadata_file = f"{metadata_type}.json"
    final_path = dst_dir / metadata_file
    is_refseq = genome_data["assembly"]["accession"].startswith("GCF_")
    # Get the sequence regions from the report and gbff files, and merge them
    report_regions = get_report_regions(report_file, is_refseq)
    if gbff_file:
        gbff_regions = get_gbff_seq_regions(gbff_file)
        seq_regions = merge_seq_regions(report_regions, gbff_regions)
    else:
        seq_regions = list(report_regions.values())
    # Exclude seq_regions from a list
    if to_exclude is not None:
        seq_regions = exclude_seq_regions(seq_regions, to_exclude)
    # Setup the BRC4_seq_region_name
    if brc4_mode:
        seq_regions = add_insdc_seq_region_name(seq_regions)
    # Add translation and mitochondrial codon tables
    add_translation_table(seq_regions)
    add_mitochondrial_codon_table(seq_regions, genome_data["species"]["taxonomy_id"])
    # Print out the file
    print_json(final_path, seq_regions)


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by the entry point of this module."""

    genome_file = argschema.fields.InputFile(
        required=True, metadata={"description": "Genome metadata JSON file path"}
    )
    report_file = argschema.fields.InputFile(
        required=True, metadata={"description": "INSDC/RefSeq sequences report file path to parse"}
    )
    gbff_file = argschema.fields.InputFile(
        required=False,
        allow_none=True,
        dump_default=None,
        metadata={"description": "INSDC/RefSeq GBFF file path to parse"},
    )
    dst_dir = argschema.fields.OutputDir(
        required=False,
        dump_default=".",
        metadata={"description": "Output folder for the processed sequence regions JSON file"},
    )
    brc4_mode = argschema.fields.Boolean(
        required=False, dump_default=True, metadata={"description": "Include INSDC sequence region names"}
    )
    to_exclude = argschema.fields.String(
        required=False,
        dump_default="None",
        metadata={"description": "Sequence region names to exclude, in a list-like string"},
    )

    @post_load
    def reformat_args(self, data: Dict[str, Any], **kwargs):
        """Processes arguments to may need additional parsing after being correctly loaded.

        Args:
            data: Loaded arguments.
            **kwargs: Arbitrary keyword arguments.

        """
        # "to_exclude" argument will be a list-like string that needs to be parsed into a list
        data["to_exclude"] = ast.literal_eval(data["to_exclude"])
        return data


def main() -> None:
    """Module's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    prepare_seq_region_metadata(
        mod.args["genome_file"],
        mod.args["report_file"],
        mod.args["gbff_file"],
        mod.args["dst_dir"],
        mod.args["brc4_mode"],
        mod.args["to_exclude"],
    )
