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
    "GFFParserError",
    "GFFSimplifier",
]

from collections import Counter
import json
import logging
from os import PathLike
from pathlib import Path
import re
from typing import Any, Dict, List, Optional

from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from importlib_resources import files

import ensembl.io.genomio.data.gff3
from ensembl.io.genomio.utils.json_utils import get_json
from .extract_annotation import FunctionalAnnotations
from .id_allocator import StableIDAllocator


class Records(list):
    """List of GFF3 SeqRecords."""

    def to_gff(self, out_gff_path: PathLike) -> None:
        """Print out the current list of records in a GFF3 file.

        Args:
            out_gff_path: Path to GFF3 file where to write the records.
        """
        with Path(out_gff_path).open("w") as out_gff_fh:
            GFF.write(self, out_gff_fh)


class GFFParserError(Exception):
    """Error when parsing a GFF3 file."""


class GFFSimplifier:
    """Parse a GGF3 file and output a cleaned up GFF3 + annotation json file.

    Raises:
        GFFParserError: If an error cannot be automatically fixed.
    """

    # pylint: disable=too-many-public-methods
    # Some rework needed at some point

    # Multiple parameters to automate various fixes
    skip_unrecognized = True
    gene_cds_skip_others = False
    allow_pseudogene_with_CDS = False
    exclude_seq_regions: List = []
    fail_types: Dict[str, int] = {}
    stable_ids = StableIDAllocator()

    def __init__(self, genome_path: Optional[PathLike] = None):
        biotypes_json = files(ensembl.io.genomio.data.gff3) / "biotypes.json"
        self._biotypes = get_json(biotypes_json)
        self.records = Records()
        self.genome: Dict[str, Dict[str, Any]] = {}
        if genome_path:
            with Path(genome_path).open("r") as genome_fh:
                self.genome = json.load(genome_fh)
        self.annotations = FunctionalAnnotations(self.genome)
        self._set_id_prefix()

    def _set_id_prefix(self) -> None:
        """Sets the ID prefix using the organism abbrev if it exists in the genome metadata."""
        try:
            org = self.genome["BRC4"]["organism_abbrev"]
        except KeyError:
            prefix = "TMP_PREFIX_"
        else:
            prefix = "TMP_" + org + "_"
        self.stable_ids.prefix = prefix

    def simpler_gff3(self, in_gff_path: PathLike) -> None:
        """Loads a GFF3 from INSDC and rewrites it in a simpler version, whilst also writing a
        functional annotation file.
        """

        with Path(in_gff_path).open("r") as in_gff_fh:
            for record in GFF.parse(in_gff_fh):
                if record.id in self.exclude_seq_regions:
                    logging.debug(f"Skip seq_region {record.id}")
                    continue

                # Clean all root features and make clean record
                clean_record = SeqRecord(record.seq, id=record.id)
                for feature in record.features:
                    clean_feature = self.simpler_gff3_feature(feature)
                    if clean_feature is not None:
                        clean_record.features.append(clean_feature)
                self.records.append(clean_record)

            if self.fail_types:
                fail_errors = "\n   ".join(self.fail_types.keys())
                if self.skip_unrecognized:
                    raise GFFParserError(f"Unrecognized types found:\n   {fail_errors}")

    def simpler_gff3_feature(self, feat: SeqFeature) -> Optional[SeqFeature]:
        """Creates a simpler version of a GFF3 feature.

        If the feature is invalid/skippable, returns None.

        """

        ignored_gene_types = self._biotypes["gene"]["ignored"]
        allowed_non_gene_types = self._biotypes["non_gene"]["supported"]
        allowed_gene_types = self._biotypes["gene"]["supported"]
        transcript_types = self._biotypes["transcript"]["supported"]

        # Skip explictly ignored features
        if feat.type in ignored_gene_types:
            return None

        # Special processing of non-gene features
        if feat.type in allowed_non_gene_types:
            if feat.type in ("mobile_genetic_element", "transposable_element"):
                feat = self.format_mobile_element(feat)
                return feat
            return None

        # From here we expect only genes
        gene = feat

        if gene.type == "protein_coding_gene":
            gene.type = "gene"

        # Create actual genes from transcripts/CDS top level features
        if gene.type in transcript_types:
            gene = self.transcript_gene(gene)
        elif gene.type == "CDS":
            gene = self.cds_gene(gene)

        # What to do with unsupported gene types
        if gene.type not in allowed_gene_types:
            self.fail_types["gene=" + gene.type] = 1
            logging.debug(f"Unsupported gene type: {gene.type} (for {gene.id})")
            if self.skip_unrecognized:
                return None

        # Normalize, store annotation, and return the cleaned up gene
        gene = self.normalize_gene(gene)
        self.store_gene_annotations(gene)
        return self.clean_gene(gene)

    def store_gene_annotations(self, gene: SeqFeature) -> None:
        """Record the functional_annotations of the gene and its children features."""
        self.annotations.add_feature(gene, "gene")

        cds_found = False
        for transcript in gene.sub_features:
            self.annotations.add_feature(transcript, "transcript", gene.id)
            for feat in transcript.sub_features:
                if feat.type != "CDS":
                    continue
                # Store CDS functional annotation only once
                if not cds_found:
                    cds_found = True
                    self.annotations.add_feature(feat, "translation", transcript.id)

    def clean_gene(self, gene: SeqFeature) -> SeqFeature:
        """Return the same gene without qualifiers unrelated to the gene structure."""

        old_gene_qualifiers = gene.qualifiers
        try:
            gene.qualifiers = {"ID": gene.id, "source": old_gene_qualifiers["source"]}
        except KeyError as err:
            raise KeyError(f"Missing source for {gene.id}") from err
        for transcript in gene.sub_features:
            # Replace qualifiers
            old_transcript_qualifiers = transcript.qualifiers
            transcript.qualifiers = {
                "ID": transcript.id,
                "Parent": gene.id,
            }
            if "source" in old_transcript_qualifiers:
                transcript.qualifiers["source"] = old_transcript_qualifiers["source"]

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

    # FORMATTERS
    def format_mobile_element(self, feat: SeqFeature) -> SeqFeature:
        """Given a mobile_genetic_element feature, transform it into a transposable_element"""

        # Change mobile_genetic_element into a transposable_element feature
        if feat.type == "mobile_genetic_element":
            mobile_element_type = feat.qualifiers.get("mobile_element_type", [])
            if mobile_element_type:
                # Get the type (and name) from the attrib
                if ":" in mobile_element_type[0]:
                    element_type, element_name = mobile_element_type[0].split(":")
                    description = f"{element_type} ({element_name})"
                else:
                    element_type = mobile_element_type[0]
                    description = element_type

                # Keep the metadata in the description if the type is known
                if element_type in ("transposon", "retrotransposon"):
                    feat.type = "transposable_element"
                    if not feat.qualifiers.get("product"):
                        feat.qualifiers["product"] = [description]
                else:
                    logging.warning(
                        f"Mobile genetic element 'mobile_element_type' is not transposon: {element_type}"
                    )
                    return feat
            else:
                logging.warning("Mobile genetic element does not have a 'mobile_element_type' tag")
                return feat
        elif feat.type == "transposable_element":
            pass
        else:
            logging.warning(f"Feature {feat.id} is not a supported TE feature {feat.type}")
            return feat

        # Generate ID if needed and add it to the functional annotation
        feat.id = self.stable_ids.normalize_gene_id(feat)
        self.annotations.add_feature(feat, "transposable_element")
        feat.qualifiers = {"ID": feat.id}

        return feat

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

    def normalize_gene(self, gene: SeqFeature) -> SeqFeature:
        """Returns a normalized gene structure, separate from the functional elements.

        Args:
            gene: Gene object to normalize.
            functional_annotation: List of feature annotations (appended by this method).

        """

        # New gene ID
        gene.id = self.stable_ids.normalize_gene_id(gene)

        # Gene with no subfeatures: need to create a transcript at least
        if len(gene.sub_features) == 0:
            logging.debug(f"Insert transcript for lone gene {gene.id}")
            transcript = self.transcript_for_gene(gene)
            gene.sub_features = [transcript]

        # Count features
        fcounter = Counter([feat.type for feat in gene.sub_features])

        # Transform gene - CDS to gene-transcript-exon-CDS
        if len(fcounter) == 1:
            if fcounter.get("CDS"):
                num_subs = len(gene.sub_features)
                logging.debug(f"Insert transcript-exon feats for {gene.id} ({num_subs} CDSs)")
                transcripts = self.gene_to_cds(gene)
                gene.sub_features = transcripts

            # Transform gene - exon to gene-transcript-exon
            elif fcounter.get("exon"):
                num_subs = len(gene.sub_features)
                logging.debug(f"Insert transcript for {gene.id} ({num_subs} exons)")
                transcript = self.gene_to_exon(gene)
                gene.sub_features = [transcript]
        else:
            # Check that we don't mix
            if fcounter.get("mRNA") and fcounter.get("CDS"):
                # Move CDS(s) from parent gene to parent mRNA if needed
                gene = self.move_cds_to_mrna(gene)
            if fcounter.get("mRNA") and fcounter.get("exon"):
                # Special case with extra exons
                gene = self.clean_extra_exons(gene)

        # Remove CDS from pseudogenes
        if gene.type == "pseudogene" and not self.allow_pseudogene_with_CDS:
            self.remove_cds_from_pseudogene(gene)

        # TRANSCRIPTS
        gene = self._normalize_transcripts(gene)

        # PSEUDOGENE CDS IDs
        if gene.type == "pseudogene" and self.allow_pseudogene_with_CDS:
            self.stable_ids.normalize_pseudogene_cds_id(gene)

        return gene

    def _normalize_transcripts(self, gene: SeqFeature) -> SeqFeature:
        """Returns a normalized transcript."""

        allowed_transcript_types = self._biotypes["transcript"]["supported"]
        ignored_transcript_types = self._biotypes["transcript"]["ignored"]
        skip_unrecognized = self.skip_unrecognized

        transcripts_to_delete = []
        for count, transcript in enumerate(gene.sub_features):
            if (
                transcript.type not in allowed_transcript_types
                and transcript.type not in ignored_transcript_types
            ):
                self.fail_types["transcript=" + transcript.type] = 1
                logging.warning(
                    f"Unrecognized transcript type: {transcript.type}" f" for {transcript.id} ({gene.id})"
                )
                if skip_unrecognized:
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

        return gene

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

                self.fail_types[f"sub_transcript={feat.type}"] = 1
                logging.warning(
                    f"Unrecognized exon type for {feat.type}: {feat.id}"
                    f" (for transcript {transcript.id} of type {transcript.type})"
                )
                if self.skip_unrecognized:
                    exons_to_delete.append(tcount)
                    continue

        if exons_to_delete:
            for elt in sorted(exons_to_delete, reverse=True):
                transcript.sub_features.pop(elt)
        return transcript

    # COMPLETION
    def transcript_gene(self, ncrna: SeqFeature) -> SeqFeature:
        """Create a gene for lone transcripts: 'gene' for tRNA/rRNA, and 'ncRNA' for all others

        Args:
            ncrna: the transcript for which we want to create a gene.

        Returns:
            The gene that contains the transcript.

        """
        new_type = "ncRNA_gene"
        if ncrna.type in ("tRNA", "rRNA"):
            new_type = "gene"
        logging.debug(f"Put the transcript {ncrna.type} in a {new_type} parent feature")
        gene = SeqFeature(ncrna.location, type=new_type)
        gene.qualifiers["source"] = ncrna.qualifiers["source"]
        gene.sub_features = [ncrna]
        gene.id = ncrna.id

        return gene

    def cds_gene(self, cds: SeqFeature) -> SeqFeature:
        """Returns a gene created for a lone CDS."""

        logging.debug(f"Put the lone CDS in gene-mRNA parent features for {cds.id}")

        # Create a transcript, add the CDS
        transcript = SeqFeature(cds.location, type="mRNA")
        transcript.qualifiers["source"] = cds.qualifiers["source"]
        transcript.sub_features = [cds]

        # Add an exon too
        exon = SeqFeature(cds.location, type="exon")
        exon.qualifiers["source"] = cds.qualifiers["source"]
        transcript.sub_features.append(exon)

        # Create a gene, add the transcript
        gene_type = "gene"
        if ("pseudo" in cds.qualifiers) and (cds.qualifiers["pseudo"][0] == "true"):
            gene_type = "pseudogene"
        gene = SeqFeature(cds.location, type=gene_type)
        gene.qualifiers["source"] = cds.qualifiers["source"]
        gene.sub_features = [transcript]
        gene.id = self.stable_ids.generate_gene_id()

        return gene

    def transcript_for_gene(self, gene: SeqFeature) -> SeqFeature:
        """Returns a transcript, from a gene without one."""

        transcript = SeqFeature(gene.location, type="transcript")
        transcript.qualifiers["source"] = gene.qualifiers["source"]
        transcript.sub_features = []

        return transcript

    def gene_to_cds(self, gene: SeqFeature) -> List[SeqFeature]:
        """Returns a list of transcripts (with exons), from a gene with only CDS children."""

        gene_cds_skip_others = self.gene_cds_skip_others
        transcripts_dict = {}
        del_transcript = []

        for count, cds in enumerate(gene.sub_features):
            if cds.type != "CDS":
                if gene_cds_skip_others:
                    del_transcript.append(count)
                    continue
                raise GFFParserError(
                    "Can not create a chain 'transcript - exon - CDS'"
                    f" when the gene children are not all CDSs"
                    f" ({cds.id} of type {cds.type} is child of gene {gene.id})"
                )

            exon = SeqFeature(cds.location, type="exon")

            # Add to transcript or create a new one
            if cds.id not in transcripts_dict:
                logging.debug(f"Create new mRNA for {cds.id}")
                transcript = self.build_transcript(gene)
                transcripts_dict[cds.id] = transcript
            exon.qualifiers["source"] = gene.qualifiers["source"]
            transcripts_dict[cds.id].sub_features.append(exon)
            transcripts_dict[cds.id].sub_features.append(cds)

        for elt in sorted(del_transcript, reverse=True):
            gene.sub_features.pop(elt)

        transcripts = list(transcripts_dict.values())

        return transcripts

    def build_transcript(self, gene: SeqFeature) -> SeqFeature:
        """Returns a transcript with same metadata as the gene provided."""

        transcript = SeqFeature(gene.location, type="mRNA")
        transcript.qualifiers["source"] = gene.qualifiers["source"]
        transcript.sub_features = []
        return transcript

    def move_cds_to_mrna(self, gene: SeqFeature) -> SeqFeature:
        """Move CDS child features of a gene, to the mRNA.

        This is to fix the case where we have the following structure:
        gene -> [ mRNA, CDSs ]
        and change it to
        gene -> [ mRNA -> [ CDSs ] ]
        The mRNA might have exons, in which case check that they match the CDS coordinates.

        Raises an exception if the feature structure is not recognized.

        Args:
            A gene with only one transcript, to check and fix.

        Returns:
            The gene where the CDSs have been moved, if needed.

        """
        # First, count the types
        mrnas = []
        cdss = []

        gene_subf_clean = []
        for subf in gene.sub_features:
            if subf.type == "mRNA":
                mrnas.append(subf)
            elif subf.type == "CDS":
                cdss.append(subf)
            else:
                gene_subf_clean.append(subf)

        if len(cdss) == 0:
            # Nothing to fix here, no CDSs to move
            return gene
        if len(mrnas) > 1:
            raise GFFParserError(
                f"Can't fix gene {gene.id}: contains several mRNAs and CDSs, all children of the gene"
            )

        mrna = mrnas[0]

        # Check if there are exons (or CDSs) under the mRNA
        sub_exons = []
        sub_cdss = []
        for subf in mrna.sub_features:
            if subf.type == "CDS":
                sub_cdss.append(subf)
            elif subf.type == "exon":
                sub_exons.append(subf)

        self._check_sub_cdss(gene, sub_cdss)
        self._check_sub_exons(gene, cdss, sub_exons)

        logging.debug(f"Gene {gene.id}: move {len(cdss)} CDSs to the mRNA")
        # No more issues? move the CDSs
        mrna.sub_features += cdss
        # And remove them from the gene
        gene.sub_features = gene_subf_clean
        gene.sub_features.append(mrna)

        return gene

    @staticmethod
    def _check_sub_cdss(gene: SeqFeature, sub_cdss: List[SeqFeature]) -> None:
        if len(sub_cdss) > 0:
            raise GFFParserError(f"Gene {gene.id} has CDSs as children of the gene and mRNA")

    @staticmethod
    def _check_sub_exons(gene: SeqFeature, cdss: SeqFeature, sub_exons: List[SeqFeature]) -> None:
        """Check that the exons of the mRNA and the CDSs match"""

        if len(sub_exons) > 0:
            # Check that they match the CDS outside
            if len(sub_exons) == len(cdss):
                # Now that all coordinates are the same
                coord_exons = [f"{exon.location}" for exon in sub_exons]
                coord_cdss = [f"{cds.location}" for cds in cdss]

                if coord_exons != coord_cdss:
                    raise GFFParserError(f"Gene {gene.id} CDSs and exons under the mRNA do not match")
            else:
                raise GFFParserError(
                    f"Gene {gene.id} CDSs and exons under the mRNA do not match (different count)"
                )

    def clean_extra_exons(self, gene: SeqFeature) -> SeqFeature:
        """Remove extra exons, already existing in the mRNA.

        This is a special case where a gene contains proper mRNAs, etc. but also
        extra exons for the same features. Those exons usually have an ID starting with
        "id-", so that's what we use to detect them.
        """
        exons = []
        mrnas = []
        others = []
        for subf in gene.sub_features:
            if subf.type == "exon":
                exons.append(subf)
            elif subf.type == "mRNA":
                mrnas.append(subf)
            else:
                others.append(subf)

        if exons and mrnas:
            exon_has_id = 0
            # Check if the exon ids start with "id-", which is an indication that they do not belong here
            for exon in exons:
                if exon.id.startswith("id-"):
                    exon_has_id += 1
            if exon_has_id:
                if exon_has_id == len(exons):
                    logging.debug(f"Remove {exon_has_id} extra exons from {gene.id}")
                    gene.sub_features = mrnas
                    gene.sub_features += others
                else:
                    raise GFFParserError(f"Can't remove extra exons for {gene.id}, not all start with 'id-'")

        return gene

    def gene_to_exon(self, gene: SeqFeature) -> SeqFeature:
        """Returns an intermediary transcript for a gene with direct exon children."""

        transcript = SeqFeature(gene.location, type="mRNA")
        transcript.qualifiers["source"] = gene.qualifiers["source"]
        transcript.sub_features = []

        for exon in gene.sub_features:
            transcript.sub_features.append(exon)

        return transcript

    def remove_cds_from_pseudogene(self, gene: SeqFeature) -> None:
        """Removes the CDS from a pseudogene.

        This assumes the CDSs are sub features of the transcript or the gene.

        """
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
