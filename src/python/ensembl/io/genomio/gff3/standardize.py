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
    "transcript_for_gene",
    "build_transcript",
    "gene_to_cds",
]

from collections import Counter
import logging
from typing import List

from Bio.SeqFeature import SeqFeature

from .exceptions import GFFParserError


def standardize_gene(gene: SeqFeature) -> None:
    """Standardize the structure of a gene model."""

    transcript_for_gene(gene)
    # Count features
    fcounter = Counter([feat.type for feat in gene.sub_features])

    # Transform gene - CDS to gene-transcript-exon-CDS
    if len(fcounter) == 1:
        if fcounter.get("CDS"):
            num_subs = len(gene.sub_features)
            logging.debug(f"Insert transcript-exon feats for {gene.id} ({num_subs} CDSs)")
            gene_to_cds(gene)

        # Transform gene - exon to gene-transcript-exon
        elif fcounter.get("exon"):
            num_subs = len(gene.sub_features)
            logging.debug(f"Insert transcript for {gene.id} ({num_subs} exons)")
            transcript = gene_to_exon(gene)
            gene.sub_features = [transcript]
    else:
        # Check that we don't mix
        if fcounter.get("mRNA") and fcounter.get("CDS"):
            # Move CDS(s) from parent gene to parent mRNA if needed
            move_cds_to_mrna(gene)
        if fcounter.get("mRNA") and fcounter.get("exon"):
            # Special case with extra exons
            clean_extra_exons(gene)


def transcript_for_gene(gene: SeqFeature) -> SeqFeature:
    """Add an mRNA to a protein coding gene if it doesn't have one."""

    if len(gene.sub_features) > 0 or gene.type != "gene":
        return gene

    logging.debug(f"Insert transcript for lone gene {gene.id}")
    transcript = SeqFeature(gene.location, type="transcript")
    transcript.qualifiers["source"] = gene.qualifiers["source"]
    transcript.sub_features = []
    gene.sub_features = [transcript]

    return gene


def build_transcript(gene: SeqFeature) -> SeqFeature:
    """Returns a transcript with same metadata as the gene provided."""

    transcript = SeqFeature(gene.location, type="mRNA")
    transcript.qualifiers["source"] = gene.qualifiers["source"]
    transcript.sub_features = []
    return transcript


def gene_to_cds(gene: SeqFeature, skip_non_cds: bool = False) -> SeqFeature:
    """Add intermediate transcripts to a gene with only CDS children."""

    transcripts_dict = {}
    del_transcript = []

    for count, cds in enumerate(gene.sub_features):
        if cds.type != "CDS":
            if skip_non_cds:
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
            transcript = build_transcript(gene)
            transcripts_dict[cds.id] = transcript
        exon.qualifiers["source"] = gene.qualifiers["source"]
        transcripts_dict[cds.id].sub_features.append(exon)
        transcripts_dict[cds.id].sub_features.append(cds)

    for elt in sorted(del_transcript, reverse=True):
        gene.sub_features.pop(elt)

    transcripts = list(transcripts_dict.values())

    gene.sub_features = transcripts
    return gene


def move_cds_to_mrna(gene: SeqFeature) -> SeqFeature:
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

    _check_sub_cdss(gene, sub_cdss)
    _check_sub_exons(gene, cdss, sub_exons)

    logging.debug(f"Gene {gene.id}: move {len(cdss)} CDSs to the mRNA")
    # No more issues? move the CDSs
    mrna.sub_features += cdss
    # And remove them from the gene
    gene.sub_features = gene_subf_clean
    gene.sub_features.append(mrna)

    return gene


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


def clean_extra_exons(gene: SeqFeature) -> SeqFeature:
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


def gene_to_exon(gene: SeqFeature) -> SeqFeature:
    """Returns an intermediary transcript for a gene with direct exon children."""

    transcript = SeqFeature(gene.location, type="mRNA")
    transcript.qualifiers["source"] = gene.qualifiers["source"]
    transcript.sub_features = []

    for exon in gene.sub_features:
        transcript.sub_features.append(exon)

    return transcript
