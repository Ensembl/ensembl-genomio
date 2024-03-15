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
from .restructure import restructure_gene
from .exceptions import GFFParserError


class Records(list):
    """List of GFF3 SeqRecords."""

    def from_gff(self, in_gff_path: PathLike, excluded: Optional[List[str]] = None) -> None:
        """Load records from a GFF3 file."""
        if excluded is None:
            excluded = []
        with Path(in_gff_path).open("r") as in_gff_fh:
            for record in GFF.parse(in_gff_fh):
                if record.id in excluded:
                    logging.debug(f"Skip seq_region {record.id}")
                    continue
                clean_record = SeqRecord(record.seq, id=record.id)
                clean_record.features = record.features
                self.append(clean_record)

    def to_gff(self, out_gff_path: PathLike) -> None:
        """Print out the current list of records in a GFF3 file."""
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
        skip_unrecognized: Optional[bool] = True,
        allow_pseudogene_with_cds: Optional[bool] = False,
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
        self.annotations = FunctionalAnnotations()

    def simpler_gff3(self, in_gff_path: PathLike) -> None:
        """Loads a GFF3 from INSDC and rewrites it in a simpler version, whilst also writing a
        functional annotation file.
        """
        self.records.from_gff(in_gff_path, self.exclude_seq_regions)
        for record in self.records:
            cleaned_features = []
            for feature in record.features:
                clean_feature = self.simpler_gff3_feature(feature)
                if clean_feature is not None:
                    cleaned_features.append(clean_feature)
            record.features = cleaned_features

        if self.fail_types:
            fail_errors = "\n   ".join(self.fail_types.keys())
            logging.warning(f"Unrecognized types found:\n   {fail_errors}")
            if not self.skip_unrecognized:
                raise GFFParserError("Unrecognized types found, abort")

    def simpler_gff3_feature(self, gene: SeqFeature) -> Optional[SeqFeature]:
        """Creates a simpler version of a GFF3 feature.

        If the feature is invalid/skippable, returns None.

        """
        # Special cases
        non_gene = self.normalize_non_gene(gene)
        if non_gene:
            return non_gene
        if gene.type in self._biotypes["gene"]["ignored"]:
            return None

        # Synonym
        if gene.type == "protein_coding_gene":
            gene.type = "gene"

        # Lone sub-gene features, create a gene
        gene = self.create_gene_for_lone_transcript(gene)
        gene = self.create_gene_for_lone_cds(gene)

        # What to do with unsupported gene types
        if gene.type not in self._biotypes["gene"]["supported"]:
            self.fail_types.add(f"gene={gene.type}")
            logging.debug(f"Unsupported gene type: {gene.type} (for {gene.id})")
            return None

        # Normalize and store
        gene = self.normalize_gene(gene)
        self.annotations.store_gene(gene)
        return self.clean_gene(gene)

    # COMPLETION
    def create_gene_for_lone_transcript(self, feat: SeqFeature) -> SeqFeature:
        """Create a gene for lone transcripts: 'gene' for tRNA/rRNA, and 'ncRNA_gene' for all others

        Args:
            ncrna: the transcript for which we want to create a gene.

        Returns:
            The gene that contains the transcript.

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
        """Returns a gene created for a lone CDS."""
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
        """Special case for non-genes (currently transposable elements only).
        Returns None if not applicable
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
        raise NotImplementedError(f"Unsupported non-gene: {feat.type} for {feat.id}")  # pragma: no cover

    def _normalize_mobile_genetic_element(self, feat: SeqFeature) -> SeqFeature:
        """Normalize a mobile element if it has a mobile_element_type field."""
        try:
            mobile_element_type = feat.qualifiers["mobile_element_type"]
        except KeyError:
            logging.warning("No 'mobile_element_type' tag found")
            return feat

        # Get the type (and name) from the attrib
        if ":" in mobile_element_type[0]:
            element_type, element_name = mobile_element_type[0].split(":")
            description = f"{element_type} ({element_name})"
        else:
            element_type = mobile_element_type[0]
            description = element_type

        # Keep the metadata in the description if the type is known
        if element_type in ("transposon", "retrotransposon"):
            if not feat.qualifiers.get("product"):
                feat.qualifiers["product"] = [description]
            return feat
        raise GFFParserError(f"'mobile_element_type' is not a transposon: {element_type}")

    # GENES
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
                "source": old_transcript_qualifiers["source"]
            }

            for feat in transcript.sub_features:
                old_qualifiers = feat.qualifiers
                feat.qualifiers = {
                    "ID": feat.id,
                    "Parent": transcript.id,
                    "source": old_qualifiers["source"],
                }
                if feat.type == "CDS":
                    try:
                        feat.qualifiers["phase"] = old_qualifiers["phase"]
                    except KeyError as err:
                        raise KeyError(
                            f"Missing phase for gene {gene.type} {gene.id}, CDS {feat.id} ({old_qualifiers})"
                        ) from err

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

    # PSEUDOGENES
    def normalize_pseudogene(self, gene: SeqFeature) -> None:
        """Normalize CDSs if allowed, otherwise remove them."""
        if gene.type != "pseudogene":
            return

        if self.allow_pseudogene_with_cds:
            self.stable_ids.normalize_pseudogene_cds_id(gene)
        else:
            self.remove_cds_from_pseudogene(gene)

    def remove_cds_from_pseudogene(self, gene: SeqFeature) -> None:
        """Removes the CDS from a pseudogene.

        This assumes the CDSs are sub features of the transcript or the gene.

        """
        if gene.type != "pseudogene":
            return

        gene_subfeats = []
        for transcript in gene.sub_features:
            if transcript.type == "CDS":
                logging.debug(f"Remove pseudo CDS {transcript.id}")
                continue
            new_subfeats = []
            for feat in transcript.sub_features:
                if feat.type == "CDS":
                    logging.debug(f"Remove pseudo CDS {feat.id}")
                    continue
                new_subfeats.append(feat)
            transcript.sub_features = new_subfeats
            gene_subfeats.append(transcript)
        gene.sub_features = gene_subfeats

    # TRANSCRIPTS
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

        """
        # Change V/C_gene_segment into a its corresponding transcript names
        if transcript.type in ("C_gene_segment", "V_gene_segment"):
            standard_name = transcript.qualifiers["standard_name"][0]
            biotype = transcript.type.replace("_segment", "")
            if re.search(r"\b(immunoglobulin|ig)\b", standard_name, flags=re.IGNORECASE):
                biotype = f"IG_{biotype}"
            elif re.search(r"\bt[- _]cell\b", standard_name, flags=re.IGNORECASE):
                biotype = f"TR_{biotype}"
            else:
                logging.warning(
                    f"Unexpected 'standard_name' content for feature {transcript.id}: {standard_name}"
                )
                return transcript
            transcript.type = biotype
        return transcript

    def _normalize_transcript_subfeatures(self, gene: SeqFeature, transcript: SeqFeature) -> SeqFeature:
        """Returns a transcript with normalized sub-features."""
        ignored_transcript_types = self._biotypes["transcript"]["ignored"]
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
                if "source" in old_exon_qualifiers:
                    feat.qualifiers["source"] = old_exon_qualifiers["source"]
            elif feat.type == "CDS":
                # New CDS ID
                feat.id = self.stable_ids.normalize_cds_id(feat.id)
                if feat.id in ("", gene.id, transcript.id):
                    feat.id = f"{transcript.id}_cds"
            else:
                if feat.type in ignored_transcript_types:
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
