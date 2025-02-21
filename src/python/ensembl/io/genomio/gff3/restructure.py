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
"""Restructure a gene model to a standard representation: `gene -> [ mRNAs -> [CDSs, exons] ]`"""

__all__ = [
    "restructure_gene",
    "add_transcript_to_naked_gene",
    "move_only_cdss_to_new_mrna",
    "move_only_exons_to_new_mrna",
    "move_cds_to_existing_mrna",
    "remove_extra_exons",
    "remove_cds_from_pseudogene",
]

from collections import Counter
import logging
from typing import List

from .exceptions import GFFParserError
from .features import GFFSeqFeature


def _get_feat_counts(gene: GFFSeqFeature) -> Counter:
    return Counter([feat.type for feat in gene.sub_features])


def restructure_gene(gene: GFFSeqFeature) -> None:
    """Standardize the structure of a gene model:
    - Add a transcript if there are no children
    - Move the CDS and exons to an mRNA if they are directly under the gene

    Args:
        gene: Gene feature to restructure.

    Raises:
        GFFParserError: If there are CDSs/exons remaining under the gene after applying the fixes.
    """
    # Skip if the children of the gene look ok
    counts = _get_feat_counts(gene)
    if (len(counts) > 0) and not counts.get("CDS") and not counts.get("exon"):
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
    if counts.get("CDS") or counts.get("exon"):
        raise GFFParserError(f"Gene {gene.id} contains direct CDSs and exons children")


def add_transcript_to_naked_gene(gene: GFFSeqFeature) -> None:
    """Add an unspecific transcript to a gene without any sub features."""

    if (len(gene.sub_features) > 0) or (gene.type != "gene"):
        return

    transcript = GFFSeqFeature(gene.location, type="transcript")
    transcript.qualifiers["source"] = gene.qualifiers["source"]
    gene.sub_features = [transcript]
    logging.debug(f"Inserted 1 transcript for a lone gene {gene.id}")


def move_only_cdss_to_new_mrna(gene: GFFSeqFeature) -> None:
    """Add intermediate mRNAs to a gene with only CDS children.
    Do nothing if some sub-features are not CDS.
    """

    counts = _get_feat_counts(gene)
    if (len(counts) != 1) or not counts.get("CDS"):
        return

    transcripts_dict = {}

    for cds in gene.sub_features:
        # We create as many transcripts as there are different CDS IDs
        if cds.id not in transcripts_dict:
            logging.debug(f"Create a new mRNA for {cds.id}")
            transcript = GFFSeqFeature(gene.location, type="mRNA")
            transcript.qualifiers["source"] = gene.qualifiers["source"]
            transcripts_dict[cds.id] = transcript

        # Add the CDS to the transcript
        transcripts_dict[cds.id].sub_features.append(cds)

        # Also add an exon in the same location
        exon = GFFSeqFeature(cds.location, type="exon")
        exon.qualifiers["source"] = gene.qualifiers["source"]
        transcripts_dict[cds.id].sub_features.append(exon)

    transcripts = list(transcripts_dict.values())
    gene.sub_features = transcripts

    logging.debug(f"Insert transcript-exon feats for {gene.id} ({len(transcripts)} CDSs)")


def move_only_exons_to_new_mrna(gene: GFFSeqFeature) -> None:
    """Add an mRNA for a gene that only has exons and move the exons under the mRNA.
    No change if the gene has other sub_features than exon.
    """

    counts = _get_feat_counts(gene)
    if (len(counts) != 1) or not counts.get("exon"):
        return

    transcript = GFFSeqFeature(gene.location, type="mRNA")
    transcript.qualifiers["source"] = gene.qualifiers["source"]
    transcript.sub_features = gene.sub_features
    gene.sub_features = [transcript]

    logging.debug(f"Insert transcript for {gene.id} ({len(gene.sub_features)} exons)")


def move_cds_to_existing_mrna(gene: GFFSeqFeature) -> None:
    """Move CDS child features of a gene to the mRNA.

    This is to fix the case where we have the following structure::
        gene -> [ mRNA, CDSs ]

    and change it to::
        gene -> [ mRNA -> [ CDSs ] ]

    The mRNA itself might have exons, in which case check that they match the CDS coordinates.

    Args:
        gene: Gene feature to update.

    Raises:
        GFFParserError: If the feature structure is not recognized.
    """
    counts = _get_feat_counts(gene)
    if not counts.get("mRNA") or not counts.get("CDS"):
        return
    if counts["mRNA"] > 1:
        raise GFFParserError(
            f"Can't fix gene {gene.id}: contains several mRNAs and CDSs, all children of the gene"
        )

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

    mrna = mrnas[0]

    # Check if there are exons (or CDSs) under the mRNA
    sub_cdss = []
    sub_exons = []
    for subf in mrna.sub_features:
        if subf.type == "CDS":
            sub_cdss.append(subf)
        elif subf.type == "exon":
            sub_exons.append(subf)

    # Check sub CDSs
    if sub_cdss:
        raise GFFParserError(f"Gene {gene.id} has CDSs as children in both the gene and mRNA")

    # If there are exons, check they overlap with the CDSs
    _check_sub_exons(mrna, cdss, sub_exons)

    # No more issues? Move the CDSs, and add any new exons
    mrna.sub_features += cdss
    # And remove them from the gene
    gene.sub_features = gene_subf_clean
    gene.sub_features.append(mrna)
    logging.debug(f"Gene {gene.id}: moved {len(cdss)} CDSs to the mRNA")


def _check_sub_exons(mrna: GFFSeqFeature, cdss: List[GFFSeqFeature], sub_exons: List[GFFSeqFeature]) -> None:
    """Check that the exons of the mRNA and the CDSs match.
    If there are no exons, create them from the CDSs.
    """

    new_sub_exons = []
    if sub_exons:
        # Check that they match the CDS outside
        if len(sub_exons) == len(cdss):
            # Now that all coordinates are the same
            coord_exons = [f"{exon.location}" for exon in sub_exons]
            coord_cdss = [f"{cds.location}" for cds in cdss]

            if coord_exons != coord_cdss:
                raise GFFParserError(f"Gene CDSs and exons under the mRNA {mrna.id} do not match")
        else:
            raise GFFParserError(
                f"Gene CDSs and exons under the mRNA {mrna.id} do not match (different count)"
            )
    else:
        # No exons in the mRNA? Create them with the CDS coordinates
        for cur_cds in cdss:
            sub_exon = GFFSeqFeature(cur_cds.location, type="exon")
            new_sub_exons.append(sub_exon)
    mrna.sub_features += new_sub_exons


def remove_extra_exons(gene: GFFSeqFeature) -> None:
    """Remove duplicated exons existing in both the gene and the mRNAs.

    This is a special case where a gene contains proper mRNAs, etc. but also extra exons for the same
    features. Those exons usually have an ID starting with "id-", so that is what we use to detect them.

    Args:
        gene: Gene feature to update.

    Raises:
        GFFParserError: If not all exons of this gene start with "id-".
    """
    counts = _get_feat_counts(gene)
    if not counts.get("mRNA") and not counts.get("exon"):
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
        if exon_has_id == len(exons):
            logging.debug(f"Remove {exon_has_id} extra exons from {gene.id}")
            gene.sub_features = mrnas
            gene.sub_features += others
        else:
            raise GFFParserError(f"Can't remove extra exons for {gene.id}, not all start with 'id-'")


def remove_cds_from_pseudogene(gene: GFFSeqFeature) -> None:
    """Removes the CDSs from a pseudogene.

    This assumes the CDSs are sub features of the transcript or the gene.

    """
    if gene.type != "pseudogene":
        return

    gene_subfeats = []
    for transcript in gene.sub_features:
        if transcript.type == "CDS":
            logging.debug(f"Remove pseudo CDS {transcript.id}")
        else:
            new_subfeats = []
            for feat in transcript.sub_features:
                if feat.type == "CDS":
                    logging.debug(f"Remove pseudo CDS {feat.id}")
                else:
                    new_subfeats.append(feat)
            transcript.sub_features = new_subfeats
            gene_subfeats.append(transcript)
    gene.sub_features = gene_subfeats
