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
"""Parse a Genbank file and creates cleaned up files from it:
- DNA fasta
- Peptide fasta
- Gene models GFF3
- seq_regions json
- genome metadata json

Raises:
    GFFPArseError: If the structure of the gb file can't be parsed.
    UnsupportedData: If some data is not as expected.

Returns:
    json_output: json file with a dict that contains all genome files created.
"""

__all__ = ["GBParseError", "UnsupportedData", "GenomeFiles", "FormattedFilesGenerator"]

from collections import Counter
import json
import logging
from os import PathLike
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from BCBio import GFF
from Bio import GenBank, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


class GBParseError(Exception):
    """Error when parsing the Genbank file."""


class UnsupportedData(Exception):
    """When an expected data is not supported by the current parser."""


class GenomeFiles(dict):
    """Store the representation of the genome files created."""

    def __init__(self, out_dir: PathLike = Path.cwd()) -> None:
        super().__init__()
        out_dir = Path(out_dir)
        self["genome"] = out_dir / "genome.json"
        self["seq_region"] = out_dir / "seq_region.json"
        self["fasta_dna"] = out_dir / "dna.fasta"
        self["fasta_pep"] = out_dir / "pep.fasta"
        self["gene_models"] = out_dir / "genes.gff"


class FormattedFilesGenerator:
    """
    Contains a parser to load data from a file, and output a set of files that follow our schema
    for input into a core database
    """

    locations = {
        "mitochondrion": "mitochondrial_chromosome",
        "apicoplast": "apicoplast_chromosome",
        "chloroplast": "chloroplast_chromosome",
        "chromoplast": "chromoplast_chromosome",
        "cyanelle": "cyanelle_chromosome",
        "leucoplast": "leucoplast_chromosome",
    }

    allowed_feat_types = [
        "gene",
        "transcript",
        "tRNA",
        "rRNA",
        "CDS",
    ]

    def __init__(self, prod_name: str, gb_file: Path, prefix: str = ""):
        self.prefix = prefix
        self.seq_records: List[SeqRecord] = []
        self.prod_name = prod_name
        self.gb_file = gb_file
        self.files = GenomeFiles()

    def set_prefix(self, prefix):
        """
        Define a prefix to add to the feature IDs
        """
        if prefix:
            self.prefix = prefix

    def set_production_name(self, prod_name):
        """
        Define a production_name for the genome
        """
        if prod_name:
            self.prod_name = prod_name

    def extract_gb(self, out_dir: Optional[PathLike]) -> Dict[str, Path]:
        """Extract data from a Genbank file and create files from it."""
        if out_dir is not None:
            self.files = GenomeFiles(out_dir)
        self.set_prefix(self.prefix)
        self.set_production_name(self.prod_name)
        self.parse_genbank(Path(self.gb_file))

        # Output the gff3 file
        return self.files

    def parse_genbank(self, gb_file):
        """
        Load a sequence from a Genbank file
        """

        organella = self._get_organella(gb_file)
        logging.debug(f"Organella loaded: {organella}")

        with open(gb_file, "r") as gbh:
            for record in SeqIO.parse(gbh, "genbank"):
                # We don't want the record description (especially for the fasta file)
                record.description = ""
                record.organelle = None
                if record.id in organella:
                    record.organelle = organella[record.id]
                self.seq_records.append(record)

            self._write_genome_json()
            self._write_genes_gff()
            self._write_seq_region_json()
            self._write_fasta_dna()

    def _get_organella(self, gb_file):
        """
        Retrieve the organelle from the genbank file, using the specific GenBank object,
        because SeqIO does not support this field
        """
        organella = {}
        with open(gb_file, "r") as gbh:
            for record in GenBank.parse(gbh):
                accession = record.version
                for q in record.features[0].qualifiers:
                    if q.key == "/organelle=":
                        organelle = q.value.replace('"', "")
                        organella[accession] = organelle
        return organella

    def _write_fasta_dna(self):
        logging.debug(f"Write {len(self.seq_records)} DNA sequences to {self.files['fasta_dna']}")
        with open(self.files["fasta_dna"], "w") as fasta_fh:
            SeqIO.write(self.seq_records, fasta_fh, "fasta")

    def _write_genes_gff(self) -> None:
        """Extract gene models from the record, and write a GFF and peptide fasta file.
        Raise GBParseError If the IDs in all the records are not unique."""
        peptides = []
        records = []
        all_ids = []

        for record in self.seq_records:
            new_record, rec_ids, rec_peptides = self._parse_record(record)
            if new_record.features:
                records.append(new_record)
            all_ids += rec_ids
            peptides += rec_peptides

        logging.debug(f"Write {len(records)} gene records to {self.files['gene_models']}")
        if records:
            with self.files["gene_models"].open("w") as gff_fh:
                GFF.write(records, gff_fh)

        logging.debug(f"Write {len(peptides)} peptide sequences to {self.files['fasta_pep']}")
        if peptides:
            with self.files["fasta_pep"].open("w") as fasta_fh:
                SeqIO.write(peptides, fasta_fh, "fasta")

        logging.debug("Check that IDs are unique")
        count = dict(Counter(all_ids))
        num_duplicates = 0
        for key in count:
            if count[key] > 1:
                num_duplicates += 1
                logging.warning(f"ID {key} is duplicated {count[key]} times")
        if num_duplicates > 0:
            raise GBParseError(f"Some {num_duplicates} IDs are duplicated")

    def _parse_record(self, record: SeqRecord) -> Tuple[SeqRecord, List[SeqRecord], List[str]]:
        all_ids: List[str] = []
        peptides: List[SeqFeature] = []
        feats: Dict[str, SeqFeature] = {}

        for feat in record.features:
            # Silently skip any unsupported feature type
            if feat.type not in self.allowed_feat_types:
                continue

            # Create a clean clone of the feature
            gff_qualifiers = feat.qualifiers
            gff_feat = SeqFeature(
                location=feat.location,
                type=feat.type,
                strand=feat.location.strand,
                qualifiers=gff_qualifiers,
            )
            # Only Genes should have a name: use either attribute gene or locus_tag
            gene_name = gff_qualifiers.get("gene", [None])[0]
            if gene_name is None:
                gene_name = gff_qualifiers.get("locus_tag", [None])[0]

            # Parse this gene
            if gene_name is not None:
                gene_feats, gene_ids, gene_peptides = self._parse_gene_feat(gff_feat, gene_name)
                peptides += gene_peptides
                feats = {**feats, **gene_feats}
                all_ids += gene_ids

            # No gene ID: parse if it is a tRNA or rRNA
            elif gff_feat.type in ("tRNA", "rRNA"):
                rna_feats, rna_ids = self._parse_rna_feat(gff_feat)
                feats = {**feats, **rna_feats}
                all_ids += rna_ids

            # Any other case? Fail here and check if we shoud support it, or add it to unsupported list
            else:
                raise GBParseError(f"No ID for allowed feature: {feat}")

        new_record = SeqRecord(record.seq, record.id)
        new_record.features = feats.values()
        return new_record, all_ids, peptides

    def _parse_gene_feat(
        self, gene_feat: SeqFeature, gene_name: str
    ) -> Tuple[Dict[str, SeqFeature], List[str], List[SeqFeature]]:
        gene_id = self.prefix + gene_name
        gene_qualifiers = gene_feat.qualifiers
        new_feats: Dict[str, Any] = {}
        peptides: List[SeqFeature] = []
        all_ids: List[str] = []

        if gene_feat.type == "gene":
            if "pseudo" in gene_qualifiers:
                gene_feat.type = "pseudogene"
            gene_feat.qualifiers["ID"] = gene_id
            gene_feat.qualifiers["Name"] = gene_name
            if "gene" in gene_feat.qualifiers:
                del gene_feat.qualifiers["gene"]
            if "locus_tag" in gene_feat.qualifiers:
                del gene_feat.qualifiers["locus_tag"]
            new_feats[str(gene_id)] = gene_feat
            all_ids.append(str(gene_id))

        if gene_feat.type in ("tRNA", "rRNA"):
            tr_id = gene_id + "_t1"
            gene_feat.qualifiers["ID"] = tr_id
            gene_feat.qualifiers["Parent"] = gene_id
            if "gene" in gene_feat.qualifiers:
                del gene_feat.qualifiers["gene"]
            if "locus_tag" in gene_feat.qualifiers:
                del gene_feat.qualifiers["locus_tag"]
            new_feats[str(tr_id)] = gene_feat
            all_ids.append(str(tr_id))

        if gene_feat.type == "CDS":
            if "pseudo" in gene_qualifiers:
                gene_feat.type = "exon"
            cds_id = gene_id + "_p1"
            tr_id = gene_id + "_t1"
            gene_feat.qualifiers["ID"] = cds_id
            gene_feat.qualifiers["Parent"] = tr_id
            if "gene" in gene_feat.qualifiers:
                del gene_feat.qualifiers["gene"]
            if "locus_tag" in gene_feat.qualifiers:
                del gene_feat.qualifiers["locus_tag"]

            # Add fasta to pep fasta file
            if "translation" in gene_qualifiers:
                new_pep_record = SeqRecord(Seq(gene_qualifiers["translation"][0]), id=cds_id)
                peptides.append(new_pep_record)

            # Also create a parent transcript for this translation
            tr_qualifiers = {"ID": tr_id, "Name": gene_name, "Parent": gene_id}
            gff_tr = SeqFeature(
                location=gene_feat.location,
                type="mRNA",
                strand=gene_feat.location.strand,
                qualifiers=tr_qualifiers,
            )
            new_feats[str(tr_id)] = gff_tr
            new_feats[str(cds_id)] = gene_feat
            all_ids.append(str(tr_id))
            all_ids.append(str(cds_id))

        return new_feats, all_ids, peptides

    def _parse_rna_feat(self, rna_feat: SeqFeature) -> Tuple[Dict[str, SeqFeature], List[SeqFeature]]:
        new_feats: Dict[str, Any] = {}
        all_ids: List[str] = []

        gff_qualifiers = rna_feat.qualifiers
        feat_name = gff_qualifiers["product"][0]
        gene_id = self.prefix + feat_name

        parts = gene_id.split(" ")
        if len(parts) > 2:
            logging.info(f"Shortening gene_id to {parts[0]}")
            gene_id = parts[0]
        gene_id = self._uniquify_id(gene_id, all_ids)

        feat_id = gene_id + "_t1"
        rna_feat.qualifiers["ID"] = feat_id
        rna_feat.qualifiers["Name"] = feat_name
        rna_feat.qualifiers["Parent"] = gene_id

        # Also create a parent gene for this transcript
        gene_qualifiers = {
            "ID": gene_id,
            "Name": feat_name,
        }
        gff_gene = SeqFeature(
            location=rna_feat.location,
            type="gene",
            strand=rna_feat.location.strand,
            qualifiers=gene_qualifiers,
        )
        new_feats[str(gene_id)] = gff_gene
        new_feats[str(feat_id)] = rna_feat
        all_ids.append(str(gene_id))
        all_ids.append(str(feat_id))

        return new_feats, all_ids

    def _uniquify_id(self, gene_id, all_ids):
        """Ensure the gene id used is unique,
        and append a number otherwise, starting at 2
        """

        new_id = gene_id
        num = 1
        while new_id in all_ids:
            num += 1
            new_id = f"{gene_id}_{num}"
        if gene_id != new_id:
            logging.info(f"Make gene id unique: {gene_id} -> {new_id}")

        return new_id

    def _write_seq_region_json(self):
        json_array = []

        for seq in self.seq_records:
            codon_table = self._get_codon_table(seq)
            if codon_table is None:
                logging.warning(
                    (
                        "No codon table found. Make sure to change the codon table number in "
                        f"{self.files['seq_region']} manually if it is not the standard codon table."
                    )
                )

                codon_table = 1
            else:
                codon_table = int(codon_table)
            seq_obj = {
                "name": seq.id,
                "coord_system_level": "chromosome",
                "circular": (seq.annotations["topology"] == "circular"),
                "codon_table": codon_table,
                "length": len(seq.seq),
            }
            if seq.organelle:
                seq_obj["location"] = self._prepare_location(seq.organelle)
                if not codon_table:
                    logging.warning(
                        (
                            f"'{seq.organelle}' is an organelle: make sure to change the codon table number "
                            f"in {self.files['seq_region']} manually if it is not the standard codon table"
                        )
                    )

            # Additional attributes for Ensembl
            seq_obj["added_sequence"] = {
                "accession": seq.id,
                "assembly_provider": {
                    "name": "GenBank",
                    "url": "https://www.ncbi.nlm.nih.gov/genbank",
                },
            }
            if not seq_obj["added_sequence"]["assembly_provider"]["name"]:
                logging.warning(
                    (f"Please add the relevant provider name for the assembly in {self.files['seq_region']}")
                )
            if not seq_obj["added_sequence"]["assembly_provider"]["url"]:
                logging.warning(
                    (f"Please add the relevant provider URL for the assembly in {self.files['seq_region']}")
                )

            json_array.append(seq_obj)
        with open(self.files["seq_region"], "w") as seq_fh:
            seq_fh.write(json.dumps(json_array, indent=4))

    def _get_codon_table(self, seq) -> Optional[int]:
        """
        Look at the CDS features to see if they have a codon table
        """
        for feat in seq.features:
            if feat.type == "CDS":
                quals = feat.qualifiers
                if "transl_table" in quals:
                    return quals["transl_table"][0]
                return None
        return None

    def _prepare_location(self, organelle):
        """
        Given an organelle name, returns the SO term corresponding to its location
        """
        if organelle in self.locations:
            return self.locations[organelle]
        raise UnsupportedData(f"Unkown organelle: {organelle}")

    def _write_genome_json(self):
        """
        Write a draft for the genome json file
        Only the production_name is needed, but the rest of the fields need to be given
        for the validation of the json file
        """

        prod_name = self.prod_name if self.prod_name else ""

        genome_data = {
            "species": {
                "production_name": prod_name,
                "taxonomy_id": 0,
            },
            "assembly": {"accession": "GCA_000000000", "version": 1},
            "added_seq": {},
        }

        if not genome_data["species"]["production_name"]:
            logging.warning(
                f"Please add the relevant production_name for this genome in {self.files['genome']}"
            )

        ids = [seq.id for seq in self.seq_records]
        genome_data["added_seq"]["region_name"] = ids

        with open(self.files["genome"], "w") as genome_fh:
            genome_fh.write(json.dumps(genome_data, indent=4))


def main() -> None:
    """Main script entry-point."""
    parser = ArgumentParser(description="Parse a GenBank file and create cleaned up files from it.")
    parser.add_argument_src_path("--gb_file", required=True, help="sequence accession file")
    parser.add_argument("--prefix", required=True, help="prefix to add to every feature ID")
    parser.add_argument("--prod_name", required=True, help="production name for the species")
    parser.add_argument_dst_path(
        "--out_dir", default=Path.cwd(), help="output folder where the generated files will be stored"
    )
    parser.add_log_arguments(add_log_file=True)
    args = parser.parse_args()
    init_logging_with_args(args)

    gb_extractor = FormattedFilesGenerator(prefix=args.prefix, prod_name=args.prod_name, gb_file=args.gb_file)
    gb_extractor.extract_gb(args.out_dir)


if __name__ == "__main__":
    main()
