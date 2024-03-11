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
"""Standardize the gene model representation of SeqFeatures.
"""

__all__ = [
    "standardize_gene",
    "add_transcript_to_naked_gene",
    "move_only_cdss_to_new_mrna",
    "move_only_exons_to_new_mrna",
    "move_cds_to_existing_mrna",
    "remove_extra_exons",
]

from collections import Counter
import logging
from typing import List

from Bio.SeqFeature import SeqFeature

from .exceptions import GFFParserError


def _get_feat_counts(gene: SeqFeature):
    return Counter([feat.type for feat in gene.sub_features])


def standardize_gene(gene: SeqFeature) -> None:
    """Standardize the structure of a gene model:
    - Add a transcript if there are no children
    - Move the CDS and exons to an mRNA if they are directly under the gene

    Exception:
    GFFParserError if there are CDS/exons remaining under the gene after fixes
    """
    # Skip if the children of the gene look ok
    counts = _get_feat_counts(gene)
    if len(counts) > 0 and not counts.get("CDS") and not counts.get("exon"):
        return
    
    # Make sure the gene has a transcript if nothing else
    add_transcript_to_naked_gene(gene)

    # Corrections if there are CDSs or exons directly under the gene level
    move_only_cdss_to_new_mrna(gene)
    move_only_exons_to_new_mrna(gene)
    move_cds_to_existing_mrna(gene)
    remove_extra_exons(gene)

    # Check again after fixes that no CDS or exon remain under the gene
    counts = _get_feat_counts(gene)
    if counts.get("CDS") or counts("exon"):
        GFFParserError(f"Gene {gene.id} contains direct CDSs and exons children")


def add_transcript_to_naked_gene(gene: SeqFeature) -> None:
    """Add a transcript to a gene without any sub features."""

    if len(gene.sub_features) > 0 or gene.type != "gene":
        return gene

    transcript = SeqFeature(gene.location, type="transcript")
    transcript.qualifiers["source"] = gene.qualifiers["source"]
    transcript.sub_features = []
    gene.sub_features = [transcript]
    logging.debug(f"Inserted 1 transcript for a lone gene {gene.id}")


def move_only_cdss_to_new_mrna(gene: SeqFeature) -> None:
    """Add intermediate mRNAs to a gene with only CDS children.
    Do nothing if some sub-features are not CDS.
    """

    counts = _get_feat_counts(gene)
    if len(counts) != 1 or not counts.get("CDS"):
        return

    transcripts_dict = {}

    for cds in gene.sub_features:
        # We create as many transcripts as there are different CDS IDs
        if cds.id not in transcripts_dict:
            logging.debug(f"Create a new mRNA for {cds.id}")
            transcript = SeqFeature(gene.location, type="mRNA")
            transcript.qualifiers["source"] = gene.qualifiers["source"]
            transcript.sub_features = []
            transcripts_dict[cds.id] = transcript
        
        # Add the CDS to the transcript
        transcripts_dict[cds.id].sub_features.append(cds)

        # Also add an exon in the same location
        exon = SeqFeature(cds.location, type="exon")
        exon.qualifiers["source"] = gene.qualifiers["source"]
        transcripts_dict[cds.id].sub_features.append(exon)

    transcripts = list(transcripts_dict.values())
    gene.sub_features = transcripts

    logging.debug(f"Insert transcript-exon feats for {gene.id} ({len(transcripts)} CDSs)")


def move_only_exons_to_new_mrna(gene: SeqFeature) -> None:
    """Add an mRNA for a gene that only has exons and move the exons under the mRNA.
    No change if the gene has other sub_features than exon.
    """

    counts = _get_feat_counts(gene)
    if len(counts) != 1 or not counts.get("exon"):
        return

    transcript = SeqFeature(gene.location, type="mRNA")
    transcript.qualifiers["source"] = gene.qualifiers["source"]
    transcript.sub_features = []

    for feat in gene.sub_features:
        if feat.type != "exon":
            return
        transcript.sub_features.append(feat)
    gene.sub_features = [transcript]

    logging.debug(f"Insert transcript for {gene.id} ({len(gene.sub_features)} exons)")


def move_cds_to_existing_mrna(gene: SeqFeature) -> None:
    """Move CDS child features of a gene, to the mRNA.

    This is to fix the case where we have the following structure:
    gene -> [ mRNA, CDSs ]
    and change it to
    gene -> [ mRNA -> [ CDSs ] ]
    The mRNA might have exons, in which case check that they match the CDS coordinates.

    Raises an exception if the feature structure is not recognized.
    """
    counts = _get_feat_counts(gene)
    if len(counts) != 2 and not counts.get("mRNA") and not counts.get("CDS"):
        return

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

    _check_sub_cdss(gene, sub_cdss)
    _check_sub_exons(gene, cdss, sub_exons)

    logging.debug(f"Gene {gene.id}: move {len(cdss)} CDSs to the mRNA")
    # No more issues? move the CDSs
    mrna.sub_features += cdss
    # And remove them from the gene
    gene.sub_features = gene_subf_clean
    gene.sub_features.append(mrna)


def _check_sub_cdss(gene: SeqFeature, sub_cdss: List[SeqFeature]) -> None:
    if len(sub_cdss) > 0:
        raise GFFParserError(f"Gene {gene.id} has CDSs as children of the gene and mRNA")


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


def remove_extra_exons(gene: SeqFeature) -> SeqFeature:
    """Remove duplicated exons existing in both the gene and the mRNAs.

    This is a special case where a gene contains proper mRNAs, etc. but also
    extra exons for the same features. Those exons usually have an ID starting with
    "id-", so that's what we use to detect them.
    """
    counts = _get_feat_counts(gene)
    if len(counts) != 2 and not counts.get("mRNA") and not counts.get("exon"):
        return

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