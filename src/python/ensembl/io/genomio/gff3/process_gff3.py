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
in a separate file."""


from collections import Counter
from os import PathLike
from pathlib import Path
import re
from typing import Dict, List, Optional

import json
import argschema

from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from ensembl.io.genomio.gff3.functional_annotation import FunctionalAnnotations


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


class GFFParserCommon:
    """Heritable class to share the list of feature types supported or ignored by the parser"""

    # Supported gene level biotypes
    gene_types = [
        "gene",
        "pseudogene",
        "ncRNA_gene",
    ]

    # Exception: non_gene_types that are nonetheless kept in the GFF3 file
    non_gene_types = [
        "transposable_element",
    ]

    # Supported transcript level biotypes
    transcript_types = [
        "C_gene_segment",
        "guide_RNA",
        "lnc_RNA",
        "miRNA",
        "misc_RNA",
        "mRNA",
        "ncRNA",
        "piRNA",
        "pseudogenic_rRNA",
        "pseudogenic_transcript",
        "pseudogenic_tRNA",
        "ribozyme",
        "RNase_MRP_RNA",
        "RNase_P_RNA",
        "rRNA",
        "scRNA",
        "snoRNA",
        "snRNA",
        "SRP_RNA",
        "telomerase_RNA",
        "transcript",
        "tRNA",
        "V_gene_segment",
    ]

    # Biotypes that are ignored, and removed from the final GFF3 file
    ignored_gene_types = [
        "cDNA_match",
        "centromere",
        "D_loop",
        "direct_repeat",
        "dispersed_repeat",
        "gap",
        "intron",
        "inverted_repeat",
        "long_terminal_repeat",
        "microsatellite",
        "origin_of_replication",
        "region",
        "repeat_region",
        "satellite_DNA",
        "sequence_feature",
        "sequence_secondary_structure",
        "sequence_uncertainty",
        "STS",
        "tandem_repeat",
        "telomere",
        "terminal%2Cinverted",
    ]

    # Ignored biotypes that are under a gene parent feature
    ignored_transcript_types = [
        "3'UTR",
        "5'UTR",
        "antisense_RNA",
        "intron",
        "non_canonical_five_prime_splice_site",
        "non_canonical_three_prime_splice_site",
    ]


class GFFGeneMerger(GFFParserCommon):
    """Specialized class to merge split genes in a GFF3 file, prior to further parsing."""

    def merge(self, in_gff_path: PathLike, out_gff_path: PathLike) -> int:
        """
        Merge genes in a gff that are split in multiple lines
        """
        to_merge = []
        merged: List[str] = []

        with Path(in_gff_path).open("r") as in_gff_fh, Path(out_gff_path).open("w") as out_gff_fh:
            for line in in_gff_fh:
                # Skip comments
                if line.startswith("#"):
                    out_gff_fh.write(line)
                else:
                    # Parse one line
                    line = line.rstrip()
                    fields = line.split("\t")
                    attr_fields = fields[8].split(";")
                    attrs = {}
                    for a in attr_fields:
                        (key, value) = a.split("=")
                        attrs[key] = value

                    # Check this is a gene to merge; cache it then
                    if fields[2] in self.gene_types and ("part" in attrs or "is_ordered" in attrs):
                        to_merge.append(fields)

                    # If not, merge previous gene if needed, and print the line
                    else:
                        if to_merge:
                            merged_str = []
                            for line_to_merge in to_merge:
                                merged_str.append("\t".join(line_to_merge))
                            merged.append("\n".join(merged_str) + "\n")

                            new_line = self._merge_genes(to_merge)
                            out_gff_fh.write(new_line)
                            to_merge = []
                        out_gff_fh.write(line + "\n")

            # Print last merged gene if there is one
            if to_merge:
                merged_str = []
                for line_to_merge in to_merge:
                    merged_str.append("\t".join(line_to_merge))
                merged.append("\n".join(merged_str) + "\n")

                new_line = self._merge_genes(to_merge)
                out_gff_fh.write(new_line)

        return len(merged)

    def _merge_genes(self, to_merge: List) -> str:
        """Returns a single gene gff3 line merged from separate parts.

        Args:
            to_merge: List of gff3 lines with gene parts.

        """
        print(f"Merge gene in {len(to_merge)} parts")
        min_start = -1
        max_end = -1
        for gene in to_merge:
            print(f"Merge part: {gene[8]}")
            start = int(gene[3])
            end = int(gene[4])

            if start < min_start or min_start < 0:
                min_start = start
            if end > max_end or max_end < 0:
                max_end = end

        # Take the first line as template and replace things
        new_gene = to_merge[0]
        new_gene[3] = str(min_start)
        new_gene[4] = str(max_end)

        attrs = new_gene[8]
        attrs = attrs.replace(";is_ordered=true", "")
        attrs = re.sub(r";part=\d+/\d+", "", attrs)
        new_gene[8] = attrs

        return "\t".join(new_gene) + "\n"


class GFFSimplifier(GFFParserCommon):
    """Parse a GGF3 file and output a cleaned up GFF3 + annotation json file.

    Raises:
        GFFParserError: If an error cannot be automatically fixed.
    """

    # Multiple parameters to automate various fixes
    skip_unrecognized = False
    gene_cds_skip_others = False
    allow_pseudogene_with_CDS = False
    exclude_seq_regions: List = []
    validate_gene_id = True
    min_id_length = 8
    stable_id_prefix = None
    current_stable_id_number: int = 0

    def __init__(self, genome_path: Optional[PathLike] = None, make_missing_stable_ids: bool = False):
        self.records = Records()
        self.annotations = FunctionalAnnotations()
        self.genome = {}
        if genome_path:
            with Path(genome_path).open("r") as genome_fh:
                self.genome = json.load(genome_fh)
        self.make_missing_stable_ids: bool = make_missing_stable_ids

    def simpler_gff3(self, in_gff_path: PathLike) -> None:
        """
        Load a GFF3 from INSDC and rewrite it in a simpler version,
        and also write a functional_annotation file
        """

        allowed_gene_types = self.gene_types
        ignored_gene_types = self.ignored_gene_types
        transcript_types = self.transcript_types
        allowed_non_gene_types = self.non_gene_types
        skip_unrecognized = self.skip_unrecognized
        to_exclude = self.exclude_seq_regions

        with Path(in_gff_path).open("r") as in_gff_fh:
            fail_types: Dict[str, int] = {}

            for record in GFF.parse(in_gff_fh):
                new_record = SeqRecord(record.seq, id=record.id)
                if record.id in to_exclude:
                    print(f"Skip seq_region {record.id}")
                    continue

                # Root features (usually genes)
                for feat in record.features:
                    # Skip or format depending on the feature type
                    if feat.type in ignored_gene_types:
                        continue
                    if feat.type in transcript_types:
                        feat = self.transcript_gene(feat)
                    elif feat.type == "CDS":
                        feat = self.cds_gene(feat)
                    elif feat.type in ("mobile_genetic_element", "transposable_element"):
                        feat = self.format_mobile_element(feat)

                    # Normalize the gene structure
                    if feat.type in allowed_gene_types:
                        feat = self.normalize_gene(feat, fail_types)
                    elif feat.type in allowed_non_gene_types:
                        pass
                    else:
                        fail_types["gene=" + feat.type] = 1
                        message = f"Unsupported feature type: {feat.type} (for {feat.id})"
                        print(message)
                        if skip_unrecognized:
                            del feat
                            continue

                    new_record.features.append(feat)
                self.records.append(new_record)

            if fail_types and not skip_unrecognized:
                fail_errors = "\n   ".join(fail_types.keys())
                raise GFFParserError(f"Unrecognized types found:\n   {fail_errors}")

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
                    print(f"Mobile genetic element 'mobile_element_type' is not transposon: {element_type}")
                    return feat
            else:
                print("Mobile genetic element does not have a 'mobile_element_type' tag")
                return feat
        elif feat.type == "transposable_element":
            pass
        else:
            print(f"Feature {feat.id} is not a supported TE feature {feat.type}")
            return feat

        # Generate ID if needed and add it to the functional annotation
        feat.id = self.normalize_gene_id(feat)
        self.annotations.add_feature(feat, "transposable_element")
        feat.qualifiers = {"ID": feat.id}

        return feat

    def format_gene_segments(self, transcript: SeqFeature) -> SeqFeature:
        """Returns the equivalent Ensembl biotype feature for gene segment transcript features.

        Supported features: "C_gene_segment" and "V_gene_segment".

        Args:
            transcript: Gene segment transcript feature.

        """
        # Change mobile_genetic_element into a transposable_element feature
        if transcript.type in ("C_gene_segment", "V_gene_segment"):
            standard_name = transcript.qualifiers["standard_name"][0]
            # Drop "_segment" from the transcript type
            biotype = transcript.type[:-8]
            if re.search(r"\b(immunoglobulin|ig)\b", standard_name, flags=re.IGNORECASE):
                biotype = f'IG_{biotype}'
            elif re.search(r"\bt[- _]cell\b", standard_name, flags=re.IGNORECASE):
                biotype = f'TR_{biotype}'
            else:
                print(f"Unexpected 'standard_name' content for feature {transcript.id}: {standard_name}")
                return transcript
            transcript.type = biotype
        else:
            print(f"Feature {transcript.id} is not a supported gene segment feature: {transcript.type}")
        return transcript

    def normalize_gene(self, gene: SeqFeature, fail_types: Dict[str, int]) -> SeqFeature:
        """Returns a normalized gene structure, separate from the functional elements.

        Args:
            gene: Gene object to normalize.
            functional_annotation: List of feature annotations (appended by this method).
            fail_types: List of feature types that are not supported (appended by this method).

        """

        # New gene ID
        gene.id = self.normalize_gene_id(gene)

        # Gene with no subfeatures: need to create a transcript at least
        if len(gene.sub_features) == 0:
            print(f"Insert transcript for lone gene {gene.id}")
            transcript = self.transcript_for_gene(gene)
            gene.sub_features = [transcript]

        # Count features
        fcounter = Counter([feat.type for feat in gene.sub_features])

        # Transform gene - CDS to gene-transcript-exon-CDS
        if len(fcounter) == 1:
            if fcounter.get("CDS"):
                num_subs = len(gene.sub_features)
                print(f"Insert transcript-exon feats for {gene.id} ({num_subs} CDSs)")
                transcripts = self.gene_to_cds(gene)
                gene.sub_features = transcripts

            # Transform gene - exon to gene-transcript-exon
            elif fcounter.get("exon"):
                num_subs = len(gene.sub_features)
                print(f"Insert transcript for {gene.id} ({num_subs} exons)")
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
        gene = self._normalize_transcripts(gene, fail_types)

        # PSEUDOGENE CDS IDs
        if gene.type == "pseudogene" and self.allow_pseudogene_with_CDS:
            self.normalize_pseudogene_cds(gene)

        # Finally, store gene functional annotation
        self.annotations.add_feature(gene, "gene")

        # replace qualifiers
        old_gene_qualifiers = gene.qualifiers
        gene.qualifiers = {"ID": gene.id, "source": old_gene_qualifiers["source"]}

        return gene

    def _normalize_transcripts(self, gene: SeqFeature, fail_types) -> SeqFeature:
        """Returns a normalized transcript."""
        allowed_transcript_types = self.transcript_types
        skip_unrecognized = self.skip_unrecognized

        transcripts_to_delete = []
        for count, transcript in enumerate(gene.sub_features):
            if (
                transcript.type not in allowed_transcript_types
                and transcript.type not in self.ignored_transcript_types
            ):
                fail_types["transcript=" + transcript.type] = 1
                message = (
                    f"Unrecognized transcript type: {transcript.type}" f" for {transcript.id} ({gene.id})"
                )
                print(message)
                if skip_unrecognized:
                    transcripts_to_delete.append(count)
                    continue

            # New transcript ID
            transcript_number = count + 1
            transcript.id = self.normalize_transcript_id(gene.id, transcript_number)

            if transcript.type in ("C_gene_segment", "V_gene_segment"):
                transcript = self.format_gene_segments(transcript)

            # Store transcript functional annotation
            self.annotations.add_feature(transcript, "transcript", gene.id)

            # Replace qualifiers
            old_transcript_qualifiers = transcript.qualifiers
            transcript.qualifiers = {
                "ID": transcript.id,
                "Parent": gene.id,
            }
            if "source" in old_transcript_qualifiers:
                transcript.qualifiers["source"] = old_transcript_qualifiers["source"]

            # EXONS AND CDS
            transcript = self._normalize_transcript_subfeatures(gene, transcript, fail_types)

        if transcripts_to_delete:
            for elt in sorted(transcripts_to_delete, reverse=True):
                gene.sub_features.pop(elt)

        return gene

    def _normalize_transcript_subfeatures(
        self, gene: SeqFeature, transcript: SeqFeature, fail_types
    ) -> SeqFeature:
        """Returns a transcript with normalized sub-features."""
        ignored_transcript_types = self.ignored_transcript_types
        cds_found = False
        exons_to_delete = []
        for tcount, feat in enumerate(transcript.sub_features):
            if feat.type == "exon":
                # Replace qualifiers
                old_exon_qualifiers = feat.qualifiers
                feat.qualifiers = {"Parent": transcript.id}
                if "source" in old_exon_qualifiers:
                    feat.qualifiers["source"] = old_exon_qualifiers["source"]
            elif feat.type == "CDS":
                # New CDS ID
                feat.id = self.normalize_cds_id(feat.id)
                if feat.id in ("", gene.id, transcript.id):
                    feat.id = f"{transcript.id}_cds"

                # Store CDS functional annotation (only once)
                if not cds_found:
                    cds_found = True
                    self.annotations.add_feature(feat, "translation", transcript.id)

                # Replace qualifiers
                feat.qualifiers = {
                    "ID": feat.id,
                    "Parent": transcript.id,
                    "phase": feat.qualifiers["phase"],
                    "source": feat.qualifiers["source"],
                }
            else:
                if feat.type in ignored_transcript_types:
                    exons_to_delete.append(tcount)
                    continue

                fail_types[f"sub_transcript={feat.type}"] = 1
                message = (
                    f"Unrecognized exon type for {feat.type}: {feat.id}"
                    f" (for transcript {transcript.id} of type {transcript.type})"
                )
                print(message)
                if self.skip_unrecognized:
                    exons_to_delete.append(tcount)
                    continue

        if exons_to_delete:
            for elt in sorted(exons_to_delete, reverse=True):
                transcript.sub_features.pop(elt)
        return transcript

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
        print(f"Put the transcript {ncrna.type} in a {new_type} parent feature")
        gene = SeqFeature(ncrna.location, type=new_type)
        gene.qualifiers["source"] = ncrna.qualifiers["source"]
        gene.sub_features = [ncrna]
        gene.id = ncrna.id

        return gene

    def cds_gene(self, cds: SeqFeature) -> SeqFeature:
        """Returns a gene created for a lone CDS."""

        print("Put the lone CDS in gene-mRNA parent features")

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
        gene.id = self.generate_stable_id()

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
                print(f"Create new mRNA for {cds.id}")
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

        print(f"Gene {gene.id}: move {len(cdss)} CDSs to the mRNA")
        # No more issues? move the CDSs
        mrna.sub_features += cdss
        # And remove them from the gene
        gene.sub_features = gene_subf_clean
        gene.sub_features.append(mrna)

        return gene

    @staticmethod
    def _check_sub_cdss(gene, sub_cdss) -> None:
        if len(sub_cdss) > 0:
            raise GFFParserError(f"Gene {gene.id} has CDSs as children of the gene and mRNA")

    @staticmethod
    def _check_sub_exons(gene, cdss, sub_exons) -> None:
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
                    print(f"Remove {exon_has_id} extra exons from {gene.id}")
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

    def normalize_gene_id(self, gene: SeqFeature) -> str:
        """Remove any unnecessary prefixes around the gene ID.

        Generate a new stable id if it is not recognized as valid.

        Args:
            gene: Gene feature to normalize.

        Returns:
            A normalized gene id.

        """

        prefixes = ["gene-", "gene:"]
        new_gene_id = self.remove_prefixes(gene.id, prefixes)

        # In case the gene id is not valid, use the GeneID
        if not self.valid_id(new_gene_id):
            print(f"Gene id is not valid: {new_gene_id}")
            qual = gene.qualifiers
            if "Dbxref" in qual:
                for xref in qual["Dbxref"]:
                    (db, value) = xref.split(":")
                    if db == "GeneID":
                        new_gene_id = f"{db}_{value}"
                        print(f"Using GeneID {new_gene_id} for stable_id instead of {gene.id}")
                        return new_gene_id

            # Make a new stable_id
            if self.make_missing_stable_ids:
                new_id = self.generate_stable_id()
                print(f"New id: {new_gene_id} -> {new_id}")
                return new_id
            raise GFFParserError(f"Can't use invalid gene id for {gene}")

        return new_gene_id

    def generate_stable_id(self) -> str:
        """Returns a new unique gene stable_id with a prefix.

        The id is made up of a prefix and a number, which is auto incremented.
        Define the prefix with the param "stable_id_prefix",
        or use the genome organism_abbrev and prepend "TMP_" to it.

        """
        if self.stable_id_prefix:
            prefix = self.stable_id_prefix
        else:
            if self.genome:
                org = self.genome.get("BRC4", {}).get("organism_abbrev")
            if org is None:
                prefix = "TMP_PREFIX_"
            else:
                prefix = "TMP_" + org + "_"
            self.stable_id_prefix = prefix

        number = self.current_stable_id_number + 1
        new_id = f"{prefix}{number}"
        self.current_stable_id_number = number

        return new_id

    def valid_id(self, name: str) -> bool:
        """Check that the format of a stable id is valid."""

        if not self.validate_gene_id:
            return True

        min_length = self.min_id_length

        # Trna (from tRNAscan)
        if re.search(r"^Trna", name):
            print(f"Stable id is a Trna from tRNA-scan: {name}")
            return False

        # Coordinates
        if re.search(r"^.+:\d+..\d+", name):
            print(f"Stable id is a coordinate: {name}")
            return False

        # Special characters
        if re.search(r"[ |]", name):
            print(f"Stable id contains special characters: {name}")
            return False

        # Min length
        if len(name) < min_length:
            print(f"Stable id is too short (<{min_length}) {name}")
            return False

        return True

    def normalize_transcript_id(self, gene_id: str, number: int) -> str:
        """Use a gene ID and a number to make a formatted transcript ID."""

        transcript_id = f"{gene_id}_t{number}"
        return transcript_id

    def normalize_cds_id(self, cds_id: str) -> str:
        """
        Check the CDS ID is proper:
        - Remove any unnecessary prefixes around the CDS ID
        - Delete the ID if it is not proper
        """

        prefixes = ["cds-", "cds:"]
        cds_id = self.remove_prefixes(cds_id, prefixes)

        # Special case: if the ID doesn't look like one, remove it
        # It needs to be regenerated
        if not self.valid_id(cds_id):
            cds_id = ""

        return cds_id

    def normalize_pseudogene_cds(self, gene: SeqFeature) -> None:
        """Ensure CDS from a pseudogene have a proper ID
        - different from the gene
        - derived from the gene if it is not proper
        """

        for transcript in gene.sub_features:
            for feat in transcript.sub_features:
                if feat.type == "CDS":
                    feat.id = self.normalize_cds_id(feat.id)
                    if feat.id in ("", gene.id):
                        feat.id = f"{transcript.id}_cds"
                        feat.qualifiers["ID"] = feat.id

    def remove_cds_from_pseudogene(self, gene: SeqFeature) -> None:
        """Remove CDS from a pseudogene
        This assumes the CDSs are sub features of the transcript or the gene
        """

        gene_subfeats = []
        for transcript in gene.sub_features:
            if transcript.type == "CDS":
                print(f"Remove pseudo CDS {transcript.id}")
                continue
            new_subfeats = []
            for feat in transcript.sub_features:
                if feat.type == "CDS":
                    print(f"Remove pseudo CDS {feat.id}")
                    continue
                new_subfeats.append(feat)
            transcript.sub_features = new_subfeats
            gene_subfeats.append(transcript)
        gene.sub_features = gene_subfeats

    def remove_prefixes(self, identifier: str, prefixes: List[str]) -> str:
        """
        Remove prefixes from an identifier if they are found
        Return the unaltered identifier otherwise
        """
        for prefix in prefixes:
            if identifier.startswith(prefix):
                identifier = identifier[len(prefix) :]
        return identifier


class InputSchema(argschema.ArgSchema):
    """Standardize the gene model representation of a GFF3 file, and extract the functional annotation
    in a separate file. Input arguments expected by this script:
    """

    in_gff_path = argschema.fields.InputFile(required=True, metadata={"description": "Input gene.gff3 path"})
    genome_data = argschema.fields.InputFile(metadata={"description": "genome.json path"})
    make_missing_stable_ids = argschema.fields.Boolean(
        default=True, metadata={"description": "Generate and add stable IDs when missing?"}
    )
    out_gff_path = argschema.fields.OutputFile(
        default="gene_models.gff3", metadata={"description": "Output gff path"}
    )
    out_func_path = argschema.fields.OutputFile(
        default="functional_annotation.json",
        metadata={"description": "Output functional_annotation.json path"},
    )
    retain_split_genes = argschema.fields.Boolean(
        default=False, metadata={"description": "Do not merge split genes automatically"}
    )


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)

    in_gff_path = mod.args["in_gff_path"]

    # Merge multiline gene features in a separate file
    interim_gff_path = Path(f"{in_gff_path}_INTERIM_MERGE")
    merger = GFFGeneMerger()
    num_merged_genes = merger.merge(in_gff_path, interim_gff_path)

    # If there are split genes, decide to merge, or just die
    if num_merged_genes > 0:
        if mod.args["retain_split_genes"]:
            raise GFFParserError("GFF contains split genes. Fix it or remove '--retain_split_genes' flag.")
        # Use the GFF with the merged genes for the next part
        in_gff_path = interim_gff_path

    # Load gff3 data and write a simpler version that follows our specifications
    # as well as a functional_annotation json file
    gff_data = GFFSimplifier(mod.args.get("genome_data"), mod.args["make_missing_stable_ids"])
    gff_data.simpler_gff3(in_gff_path)
    gff_data.records.to_gff(mod.args["out_gff_path"])
    gff_data.annotations.to_json(mod.args["out_func_path"])


if __name__ == "__main__":
    main()
