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
"""Standardize the gene model representation of a GFF3 file, and extract the functional annotation
in a separate file.
"""

__all__ = [
    "Records",
    "GFFSimplifier",
]

import json
import logging
from os import PathLike
from pathlib import Path
import re
from typing import List, Optional, Set

from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from importlib_resources import files

import ensembl.io.genomio.data.gff3
from ensembl.io.genomio.utils.json_utils import get_json
from .extract_annotation import FunctionalAnnotations
from .id_allocator import StableIDAllocator
from .restructure import restructure_gene, remove_cds_from_pseudogene
from .exceptions import GFFParserError, IgnoredFeatureError, UnsupportedFeatureError


class Records(list):
    """List of GFF3 SeqRecords."""

    def from_gff(self, in_gff_path: PathLike, excluded: Optional[List[str]] = None) -> None:
        """Loads records from a GFF3 file.

        Args:
            in_gff_path: Input GFF3 file path.
            excluded: Record IDs to not load from the GFF3 file.
        """
        if excluded is None:
            excluded = []
        with Path(in_gff_path).open("r") as in_gff_fh:
            for record in GFF.parse(in_gff_fh):
                if record.id in excluded:
                    logging.debug(f"Skip seq_region {record.id} - in exclusion list")
                    continue
                clean_record = SeqRecord(record.seq, id=record.id)
                clean_record.features = record.features
                self.append(clean_record)

    def to_gff(self, out_gff_path: PathLike) -> None:
        """Writes the current list of records in a GFF3 file.

        Args:
            out_gff_path: Path to GFF3 file where to write the records.
        """
        with Path(out_gff_path).open("w") as out_gff_fh:
            GFF.write(self, out_gff_fh)


class GFFSimplifier:
    """Parse a GGF3 file and output a cleaned up GFF3 + annotation json file.

    Raises:
        GFFParserError: If an error cannot be automatically fixed.
    """

    def __init__(
        self,
        genome_path: Optional[PathLike] = None,
        skip_unrecognized: bool = False,
        allow_pseudogene_with_cds: bool = False,
    ):
        """Create an object that simplifies `SeqFeature` objects.

        Args:
            genome_path: Genome metadata file.
            skip_unrecognized: Do not include unknown biotypes instead of raising an exception.
            allow_pseudogene_with_cds: Keep CDSs under pseudogenes that have them. Delete them otherwise.

        Raises:
            GFFParserError: If a biotype is unknown and `skip_unrecognized` is False.
        """
        self.skip_unrecognized = skip_unrecognized
        self.allow_pseudogene_with_cds = allow_pseudogene_with_cds

        # Load biotypes
        biotypes_json = files(ensembl.io.genomio.data.gff3) / "biotypes.json"
        self._biotypes = get_json(biotypes_json)

        # Load genome metadata
        self.genome = {}
        if genome_path:
            with Path(genome_path).open("r") as genome_fh:
                self.genome = json.load(genome_fh)

        # Other preparations
        self.stable_ids = StableIDAllocator()
        self.stable_ids.set_prefix(self.genome)
        self.exclude_seq_regions: List[str] = []
        self.fail_types: Set = set()

        # Init the actual data we will store
        self.records = Records()
        self.annotations = FunctionalAnnotations(self.get_provider_name())

    def get_provider_name(self) -> str:
        provider_name = "GenBank"
        if self.genome:
            try:
                provider_name = self.genome["assembly"]["provider_name"]
            except KeyError:
                if self.genome["assembly"]["accession"].startswith("GCF"):
                    provider_name = "RefSeq"
        else:
            logging.warning(f"No Provider name, using default {provider_name}")
        return provider_name

    def simpler_gff3(self, in_gff_path: PathLike) -> None:
        """Loads a GFF3 from INSDC and rewrites it in a simpler version, whilst also writing a
        functional annotation file.
        """
        self.records.from_gff(in_gff_path, self.exclude_seq_regions)
        for record in self.records:
            cleaned_features = []
            for feature in record.features:
                split_genes = self.normalize_mirna(feature)
                if split_genes:
                    cleaned_features += split_genes
                else:
                    try:
                        clean_feature = self.simpler_gff3_feature(feature)
                        cleaned_features.append(clean_feature)
                    except (UnsupportedFeatureError, IgnoredFeatureError) as err:
                        logging.debug(err.message)
            record.features = cleaned_features

        if self.fail_types:
            fail_errors = "\n   ".join(list(self.fail_types))
            logging.warning(f"Unrecognized types found:\n   {fail_errors}")
            if not self.skip_unrecognized:
                raise GFFParserError("Unrecognized types found, abort")

    def simpler_gff3_feature(self, gene: SeqFeature) -> SeqFeature:
        """Creates a simpler version of a GFF3 feature.

        Raises:
            IgnoredFeatureError: If the feature type is ignored.
            UnsupportedFeatureError: If the feature type is not supported.
        """
        # Special cases
        non_gene = self.normalize_non_gene(gene)
        if non_gene:
            return non_gene
        if gene.type in self._biotypes["gene"]["ignored"]:
            raise IgnoredFeatureError(f"Ignored type {gene.type} for {gene.id}")

        # Synonym
        if gene.type == "protein_coding_gene":
            gene.type = "gene"

        # Lone sub-gene features, create a gene
        gene = self.create_gene_for_lone_transcript(gene)
        gene = self.create_gene_for_lone_cds(gene)

        # What to do with unsupported gene types
        if gene.type not in self._biotypes["gene"]["supported"]:
            self.fail_types.add(f"gene={gene.type}")
            raise UnsupportedFeatureError(f"Unsupported type {gene.type} for {gene.id}")

        # Normalize and store
        gene = self.normalize_gene(gene)
        self.annotations.store_gene(gene)
        return self.clean_gene(gene)

    def create_gene_for_lone_transcript(self, feat: SeqFeature) -> SeqFeature:
        """Returns a gene for lone transcripts: 'gene' for tRNA/rRNA, and 'ncRNA_gene' for all others.

        Args:
            feat: The transcript for which we want to create a gene.
        """
        transcript_types = self._biotypes["transcript"]["supported"]
        if feat.type not in transcript_types:
            return feat

        new_type = "ncRNA_gene"
        if feat.type in ("tRNA", "rRNA"):
            new_type = "gene"
        logging.debug(f"Put the transcript {feat.type} in a {new_type} parent feature")
        new_gene = SeqFeature(feat.location, type=new_type)
        new_gene.qualifiers["source"] = feat.qualifiers["source"]
        new_gene.sub_features = [feat]

        # Use the transcript ID for the gene, and generate a sub ID for the transcript
        new_gene.id = feat.id
        new_gene.qualifiers["ID"] = new_gene.id
        feat.id = self.stable_ids.generate_transcript_id(new_gene.id, 1)
        feat.qualifiers["ID"] = feat.id

        return new_gene

    def create_gene_for_lone_cds(self, feat: SeqFeature) -> SeqFeature:
        """Returns a gene created for a lone CDS.

        Args:
            feat: The CDS for which we want to create a gene.
        """
        if feat.type != "CDS":
            return feat

        logging.debug(f"Put the lone CDS in gene-mRNA parent features for {feat.id}")

        # Create a transcript, add the CDS
        transcript = SeqFeature(feat.location, type="mRNA")
        transcript.qualifiers["source"] = feat.qualifiers["source"]
        transcript.sub_features = [feat]

        # Add an exon too
        exon = SeqFeature(feat.location, type="exon")
        exon.qualifiers["source"] = feat.qualifiers["source"]
        transcript.sub_features.append(exon)

        # Create a gene, add the transcript
        gene_type = "gene"
        if ("pseudo" in feat.qualifiers) and (feat.qualifiers["pseudo"][0] == "true"):
            gene_type = "pseudogene"
            del feat.qualifiers["pseudo"]
        new_gene = SeqFeature(feat.location, type=gene_type)
        new_gene.qualifiers["source"] = feat.qualifiers["source"]
        new_gene.sub_features = [transcript]
        new_gene.id = self.stable_ids.generate_gene_id()
        new_gene.qualifiers["ID"] = new_gene.id
        transcript.id = self.stable_ids.generate_transcript_id(new_gene.id, 1)
        transcript.qualifiers["ID"] = transcript.id

        return new_gene

    def normalize_non_gene(self, feat: SeqFeature) -> Optional[SeqFeature]:
        """Returns a normalised "non-gene" or `None` if not applicable.

        Only transposable elements supported at the moment.

        Args:
            feat: Feature to normalise.

        Raises:
            NotImplementedError: If the feature is a not supported non-gene.
        """

        if feat.type not in self._biotypes["non_gene"]["supported"]:
            return None
        if feat.type in ("mobile_genetic_element", "transposable_element"):
            feat.type = "transposable_element"
            feat = self._normalize_mobile_genetic_element(feat)
            # Generate ID if needed
            feat.id = self.stable_ids.normalize_gene_id(feat)
            feat.qualifiers["ID"] = feat.id

            self.annotations.add_feature(feat, "transposable_element")
            return self.clean_gene(feat)
        # This is a failsafe in case you add supported non-genes
        raise NotImplementedError(f"Unsupported non-gene: {feat.type} for {feat.id}")

    def _normalize_mobile_genetic_element(self, feat: SeqFeature) -> SeqFeature:
        """Normalize a mobile element if it has a mobile_element_type field."""
        try:
            mobile_element_type = feat.qualifiers["mobile_element_type"]
        except KeyError:
            logging.warning("No 'mobile_element_type' tag found")
            return feat

        # Get the type (and name) from the attrib
        element_type, _, element_name = mobile_element_type[0].partition(":")
        description = element_type
        if element_name:
            description += f" ({element_name})"

        # Keep the metadata in the description if the type is known
        if element_type in ("transposon", "retrotransposon"):
            if not feat.qualifiers.get("product"):
                feat.qualifiers["product"] = [description]
            return feat
        raise GFFParserError(f"'mobile_element_type' is not a transposon: {element_type}")

    def clean_gene(self, gene: SeqFeature) -> SeqFeature:
        """Return the same gene without qualifiers unrelated to the gene structure."""

        old_gene_qualifiers = gene.qualifiers
        gene.qualifiers = {"ID": gene.id, "source": old_gene_qualifiers["source"]}
        for transcript in gene.sub_features:
            # Replace qualifiers
            old_transcript_qualifiers = transcript.qualifiers
            transcript.qualifiers = {
                "ID": transcript.id,
                "Parent": gene.id,
                "source": old_transcript_qualifiers["source"],
            }

            for feat in transcript.sub_features:
                old_qualifiers = feat.qualifiers
                feat.qualifiers = {
                    "ID": feat.id,
                    "Parent": transcript.id,
                    "source": old_qualifiers["source"],
                }
                if feat.type == "CDS":
                    feat.qualifiers["phase"] = old_qualifiers["phase"]

        return gene

    def normalize_gene(self, gene: SeqFeature) -> SeqFeature:
        """Returns a normalized gene structure, separate from the functional elements.

        Args:
            gene: Gene object to normalize.
            functional_annotation: List of feature annotations (appended by this method).

        """

        gene.id = self.stable_ids.normalize_gene_id(gene)
        restructure_gene(gene)
        self.normalize_transcripts(gene)
        self.normalize_pseudogene(gene)

        return gene

    def normalize_pseudogene(self, gene: SeqFeature) -> None:
        """Normalize CDSs if allowed, otherwise remove them."""
        if gene.type != "pseudogene":
            return

        if self.allow_pseudogene_with_cds:
            self.stable_ids.normalize_pseudogene_cds_id(gene)
        else:
            remove_cds_from_pseudogene(gene)

    def normalize_transcripts(self, gene: SeqFeature) -> None:
        """Normalizes a transcript."""

        allowed_transcript_types = self._biotypes["transcript"]["supported"]
        ignored_transcript_types = self._biotypes["transcript"]["ignored"]

        transcripts_to_delete = []
        for count, transcript in enumerate(gene.sub_features):
            if (
                transcript.type not in allowed_transcript_types
                and transcript.type not in ignored_transcript_types
            ):
                self.fail_types.add(f"transcript={transcript.type}")
                logging.warning(
                    f"Unrecognized transcript type: {transcript.type}" f" for {transcript.id} ({gene.id})"
                )
                transcripts_to_delete.append(count)
                continue

            # New transcript ID
            transcript_number = count + 1
            transcript.id = self.stable_ids.generate_transcript_id(gene.id, transcript_number)

            transcript = self.format_gene_segments(transcript)

            # EXONS AND CDS
            transcript = self._normalize_transcript_subfeatures(gene, transcript)

        if transcripts_to_delete:
            for elt in sorted(transcripts_to_delete, reverse=True):
                gene.sub_features.pop(elt)

    def format_gene_segments(self, transcript: SeqFeature) -> SeqFeature:
        """Returns the equivalent Ensembl biotype feature for gene segment transcript features.

        Supported features: "C_gene_segment" and "V_gene_segment".

        Args:
            transcript: Gene segment transcript feature.

        Raises:
            GFFParserError: Missing or unexpected transcript's standard name.
        """
        if transcript.type not in ("C_gene_segment", "V_gene_segment"):
            return transcript

        # Change V/C_gene_segment into a its corresponding transcript names
        try:
            standard_name = transcript.qualifiers["standard_name"][0]
        except KeyError as err:
            raise GFFParserError(f"No standard_name for {transcript.type}") from err
        biotype = transcript.type.replace("_segment", "")
        if re.search(r"\b(immunoglobulin|ig)\b", standard_name, flags=re.IGNORECASE):
            transcript.type = f"IG_{biotype}"
        elif re.search(r"\bt[- _]cell\b", standard_name, flags=re.IGNORECASE):
            transcript.type = f"TR_{biotype}"
        else:
            raise GFFParserError(f"Unexpected 'standard_name' for {transcript.id}: {standard_name}")
        return transcript

    def _normalize_transcript_subfeatures(self, gene: SeqFeature, transcript: SeqFeature) -> SeqFeature:
        """Returns a transcript with normalized sub-features."""
        exons_to_delete = []
        exon_number = 1
        for tcount, feat in enumerate(transcript.sub_features):
            if feat.type == "exon":
                # New exon ID
                feat.id = f"{transcript.id}-E{exon_number}"
                exon_number += 1
                # Replace qualifiers
                old_exon_qualifiers = feat.qualifiers
                feat.qualifiers = {"Parent": transcript.id}
                feat.qualifiers["source"] = old_exon_qualifiers["source"]
            elif feat.type == "CDS":
                # New CDS ID
                feat.id = self.stable_ids.normalize_cds_id(feat.id)
                if feat.id in ("", gene.id, transcript.id):
                    feat.id = f"{transcript.id}_cds"
            else:
                if feat.type in self._biotypes["transcript"]["ignored"]:
                    exons_to_delete.append(tcount)
                    continue

                self.fail_types.add(f"sub_transcript={feat.type}")
                logging.warning(
                    f"Unrecognized exon type for {feat.type}: {feat.id}"
                    f" (for transcript {transcript.id} of type {transcript.type})"
                )
                exons_to_delete.append(tcount)
                continue

        if exons_to_delete:
            for elt in sorted(exons_to_delete, reverse=True):
                transcript.sub_features.pop(elt)
        return transcript

    def normalize_mirna(self, gene: SeqFeature) -> List[SeqFeature]:
        """Returns gene representations from a miRNA gene that can be loaded in an Ensembl database.

        Change the representation from the form `gene[ primary_transcript[ exon, miRNA[ exon ] ] ]`
        to `gene[ primary_transcript[ exon ] ]` and `gene[ miRNA[ exon ] ]`

        Raises:
            GFFParserError: If gene has more than 1 transcript, the transcript was not formatted
                correctly or there are unknown sub-features.
        """
        base_id = gene.id
        transcripts = gene.sub_features

        # Insert main gene first if needed
        old_gene = gene
        if gene.type == "primary_transcript":
            primary = old_gene
            gene = SeqFeature(primary.location, type="gene")
            gene.sub_features = [primary]
            gene.qualifiers = primary.qualifiers
            transcripts = gene.sub_features
            gene.id = f"{base_id}_0"
            gene.qualifiers["ID"] = gene.id

        if (len(transcripts) == 0) or (transcripts[0].type != "primary_transcript"):
            return []
        if len(transcripts) > 1:
            raise GFFParserError(f"Gene has too many sub_features for miRNA {gene.id}")

        # Passed the checks
        primary = transcripts[0]

        logging.debug(f"Formatting miRNA gene {gene.id}")

        new_genes = []
        new_primary_subfeatures = []
        num = 1
        for sub in primary.sub_features:
            if sub.type == "exon":
                new_primary_subfeatures.append(sub)
            elif sub.type == "miRNA":
                new_gene_id = f"{base_id}_{num}"
                num += 1
                new_gene = SeqFeature(sub.location, "gene", id=new_gene_id)
                new_gene.qualifiers = {"source": sub.qualifiers["source"], "ID": new_gene_id}
                new_gene.sub_features = [sub]
                new_genes.append(new_gene)
            else:
                raise GFFParserError(f"Unknown subtypes for miRNA features: {sub.id}")
        primary.sub_features = new_primary_subfeatures

        if not new_genes:
            logging.debug(f"Primary_transcript without miRNA in {gene.id}")
            all_genes = [gene]
        else:
            all_genes = [gene] + new_genes

        # Normalize like other genes
        all_genes_cleaned = []
        for new_gene in all_genes:
            new_gene = self.normalize_gene(new_gene)
            self.annotations.store_gene(new_gene)
            all_genes_cleaned.append(self.clean_gene(new_gene))
        return all_genes_cleaned
