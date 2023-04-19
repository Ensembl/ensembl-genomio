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


from collections import Counter
from pathlib import Path
import re
import tempfile
import json
import argschema

from os import PathLike
from typing import List, TextIO

# import eHive
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from ensembl.brc4.runnable.utils import print_json


class InputSchema(argschema.ArgSchema):
    """Input arguments expected by this script."""

    in_gff_path = argschema.fields.InputFile(
        metadata={"required": True, "description": "gene.gff file required"}
    )
    genome_data = argschema.fields.InputFile(
        metadata={"required": True, "description": "genome.json file required"}
    )


class process_gff3:
    gene_types = [
        "gene",
        "pseudogene",
        "ncRNA_gene",
    ]
    non_gene_types = [
        "transposable_element",
    ]
    transcript_types = [
        "transcript",
        "mRNA",
        "pseudogenic_transcript",
        "tRNA",
        "pseudogenic_tRNA",
        "pseudogenic_rRNA",
        "telomerase_RNA",
        "RNase_P_RNA",
        "SRP_RNA",
        "rRNA",
        "lnc_RNA",
        "snoRNA",
        "snRNA",
        "ncRNA",
        "miRNA",
        "ribozyme",
        "piRNA",
        "RNase_MRP_RNA",
        "misc_RNA",
    ]
    ignored_gene_types = [
        "intron",
        "region",
        "gap",
        "sequence_feature",
        "sequence_uncertainty",
        "microsatellite",
        "satellite_DNA",
        "cDNA_match",
        "STS",
        "telomere",
        "centromere",
        "repeat_region",
        "inverted_repeat",
        "tandem_repeat",
        "long_terminal_repeat",
        "dispersed_repeat",
        "terminal%2Cinverted",
    ]
    ignored_transcript_types = [
        "antisense_RNA",
        "3'UTR",
        "5'UTR",
        "intron",
    ]

    skip_unrecognized = False
    gene_cds_skip_others = False
    allow_pseudogene_with_CDS = False
    merge_split_genes = False
    exclude_seq_regions = []
    validate_gene_id = True
    min_id_length = 8
    make_missing_stable_id = False
    stable_id_prefix = None
    current_stable_id_number = None

    def merge_genes_gff(self, in_gff_path: Path, out_gff_fh: TextIO) -> None:
        """
        Merge genes in a gff that are split in multiple lines
        """
        tomerge = []
        merged = []

        with in_gff_path.open("r") as gff3_in:
            for line in gff3_in:
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
                        tomerge.append(fields)

                    # If not, merge previous gene if needed, and print the line
                    else:
                        if tomerge:
                            merged_str = []
                            for line_tomerge in tomerge:
                                merged_str.append("\t".join(line_tomerge))
                            merged.append("\n".join(merged_str) + "\n")

                            new_line = self.merge_genes(tomerge)
                            out_gff_fh.write(new_line)
                            tomerge = []
                        out_gff_fh.write(line + "\n")

            # Print last merged gene if there is one
            if tomerge:
                merged_str = []
                for line_tomerge in tomerge:
                    merged_str.append("\t".join(line_tomerge))
                merged.append("\n".join(merged_str) + "\n")

                new_line = self.merge_genes(tomerge)
                out_gff_fh.write(new_line)

        if merged and not self.merge_split_genes:
            count = len(merged)
            raise Exception("%s merged genes:\n%s\n" % (count, "\n".join(merged)))

    def merge_genes(self, tomerge: List) -> str:
        """Returns a single gene gff3 line merged from separate parts.

        Args:
            tomerge: List of gff3 lines with gene parts.

        """
        print(f"Merge gene in {len(tomerge)} parts")
        min_start = -1
        max_end = -1
        for gene in tomerge:
            print("Merge part: %s" % gene[8])
            start = int(gene[3])
            end = int(gene[4])

            if start < min_start or min_start < 0:
                min_start = start
            if end > max_end or max_end < 0:
                max_end = end

        # Take the first line as template and replace things
        new_gene = tomerge[0]
        new_gene[3] = str(min_start)
        new_gene[4] = str(max_end)

        attrs = new_gene[8]
        attrs = attrs.replace(";is_ordered=true", "")
        attrs = re.sub(r";part=\d+/\d+", "", attrs)
        new_gene[8] = attrs

        return "\t".join(new_gene) + "\n"

    def simpler_gff3(self, gff3_in: TextIO, out_gff_path: Path, out_funcann_path: Path) -> None:
        """
        Load a GFF3 from INSDC, and rewrite it in a simpler version,
        and also write a functional_annotation file
        """

        allowed_gene_types = self.gene_types
        ignored_gene_types = self.ignored_gene_types
        transcript_types = self.transcript_types
        allowed_non_gene_types = self.non_gene_types
        skip_unrecognized = self.skip_unrecognized
        to_exclude = self.exclude_seq_regions

        functional_annotation = []

        with out_gff_path.open("w") as gff3_out:
            new_records = []
            fail_types = {}

            for record in GFF.parse(gff3_in):
                new_record = SeqRecord(record.seq, id=record.id)
                if record.id in to_exclude:
                    print(f"Skip seq_region {record.id}")
                    continue

                # Root features (usually genes)
                for feat in record.features:
                    # Skip or format depending on the feature type
                    if feat.type in ignored_gene_types:
                        continue
                    elif feat.type in transcript_types:
                        feat = self.transcript_gene(feat)
                    elif feat.type == "CDS":
                        feat = self.cds_gene(feat)
                    elif feat.type in ("mobile_genetic_element", "transposable_element"):
                        feat = self.format_mobile_element(feat, functional_annotation)

                    # Normalize the gene structure
                    if feat.type in allowed_gene_types:
                        feat = self.normalize_gene(feat, functional_annotation, fail_types)
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
                new_records.append(new_record)

            if fail_types and not skip_unrecognized:
                raise Exception(f"Unrecognized types found ({' '.join(fail_types.keys())}): fail")

            GFF.write(new_records, gff3_out)

        # Write functional annotation
        functional_annotation = self.clean_functional_annotations(functional_annotation)
        print_json(out_funcann_path, functional_annotation)

    def format_mobile_element(self, feat, functional_annotation):
        """Given a mobile_genetic_element feature, transform it into a transposable_element"""
        quals = feat.qualifiers

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

                # Keep the metdata in the description if the type is known
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
        self.add_funcann_feature(functional_annotation, feat, "transposable_element")
        feat.qualifiers = {"ID": feat.id}

        return feat

    def normalize_gene(self, gene: SeqFeature, functional_annotation: List, fail_types: List) -> SeqFeature:
        """Returns a normalized gene structure, separate from the functional elements.

        Args:
            gene: Gene object to normalize.
            functional_annotation: List of feature annotations (appended by this method).
            fail_types: List of feature types that are not supported (appended by this method).

        """

        allowed_transcript_types = self.transcript_types
        ignored_transcript_types = self.ignored_transcript_types
        skip_unrecognized = self.skip_unrecognized

        # New gene ID
        gene.id = self.normalize_gene_id(gene)
        self.transfer_description(gene)

        # Gene with no subfeatures: need to create a transcript at least
        if len(gene.sub_features) == 0:
            print("Insert transcript for lone gene %s" % (gene.id))
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

            # Store transcript functional annotation
            self.add_funcann_feature(functional_annotation, transcript, "transcript")

            # Replace qualifiers
            old_transcript_qualifiers = transcript.qualifiers
            transcript.qualifiers = {
                "ID": transcript.id,
                "Parent": gene.id,
            }
            if "source" in old_transcript_qualifiers:
                transcript.qualifiers["source"] = old_transcript_qualifiers["source"]

            # EXONS AND CDS
            cds_found = False
            exons_to_delete = []
            for tcount, feat in enumerate(transcript.sub_features):
                if feat.type == "exon":
                    # Replace qualifiers
                    old_exon_qualifiers = feat.qualifiers
                    feat.qualifiers = {
                        "Parent": transcript.id,
                    }
                    if "source" in old_exon_qualifiers:
                        feat.qualifiers["source"] = old_exon_qualifiers["source"]
                elif feat.type == "CDS":
                    # New CDS ID
                    feat.id = self.normalize_cds_id(feat.id)
                    if feat.id == "" or feat.id == gene.id or feat.id == transcript.id:
                        feat.id = "%s_cds" % transcript.id

                    # Store CDS functional annotation (only once)
                    if not cds_found:
                        cds_found = True
                        self.add_funcann_feature(functional_annotation, feat, "translation")

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
                    else:
                        fail_types["sub_transcript=" + feat.type] = 1
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

        if transcripts_to_delete:
            for elt in sorted(transcripts_to_delete, reverse=True):
                gene.sub_features.pop(elt)

        # PSEUDOGENE CDS IDs
        if gene.type == "pseudogene" and self.allow_pseudogene_with_CDS:
            self.normalize_pseudogene_cds(gene)

        # Finally, store gene functional annotation
        self.add_funcann_feature(functional_annotation, gene, "gene")

        # replace qualifiers
        old_gene_qualifiers = gene.qualifiers
        gene.qualifiers = {"ID": gene.id, "source": old_gene_qualifiers["source"]}

        return gene

    def clean_functional_annotations(self, functional_annotation: List) -> List:
        """Returns the functional annotations list without putative product descriptions."""
        for feat in functional_annotation:
            if "description" in feat and not self.check_product(feat["description"]):
                del feat["description"]
        return functional_annotation

    def check_product(self, product: str) -> bool:
        """Check a product string.

        Args:
            product: a product name.

        Returns:
            True only if the string is valid, False otherwise.

        """
        excluded_names = re.compile(
            r"^(uncharacterized|putative|hypothetical|predicted)"
            r"( uncharacterized)?"
            r" protein"
            r"( of unknown function)?"
            r"( \(fragment\))?$"
        )

        if excluded_names.match(product.lower()):
            return False
        return True

    def transfer_description(self, gene: SeqFeature) -> None:
        """Transfer descriptions from transcripts/translations to gene/transcripts.

        If a gene has no description but the transcript has one, transfer it to the gene.
        Otherwise, if the translation has a description, transfer it to the transcript and gene.

        Args:
            gene: Gene to transfer descriptions.

        """
        allowed_transcript_types = self.transcript_types

        if "product" not in gene.qualifiers:
            for tran in gene.sub_features:
                if tran.type in allowed_transcript_types:
                    if "product" in tran.qualifiers:
                        description = tran.qualifiers["product"][0]
                        gene.qualifiers["product"] = [description]
                        return

                    # No transcript product, but a CDS product? Copy it to both transcript and gene
                    else:
                        for cds in tran.sub_features:
                            if cds.type == "CDS" and "product" in cds.qualifiers:
                                description = cds.qualifiers["product"][0]
                                tran.qualifiers["product"] = [description]
                                gene.qualifiers["product"] = [description]
                        # Continue transfering the translation products to the transcripts

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
        if "pseudo" in cds.qualifiers and cds.qualifiers["pseudo"][0] == "true":
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
                else:
                    raise Exception(
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
        gene -> [ mNRA -> [ CDSs ] ]
        The mRNA might have exons, in which case check that they match the CDS coordinates.

        Raises an exception if the feature structure is not recognized.

        Args:
            A gene with only one transcript, to check and fix.
        Returns:
            The gene where the CDSs have been moved, if needed.

        """
        # First, count the types
        feats = []
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
            raise Exception(
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

        if len(sub_cdss) > 0:
            raise Exception(f"Gene {gene.id} has CDSs as children of the gene and mRNA")
        if len(sub_exons) > 0:
            # Check that they match the CDS outside
            if len(sub_exons) == len(cdss):
                # Now that all coordinates are the same
                coord_exons = [f"{exon.location}" for exon in sub_exons]
                coord_cdss = [f"{cds.location}" for cds in cdss]

                if coord_exons != coord_cdss:
                    raise Exception(f"Gene {gene.id} CDSs and exons under the mRNA do not match")
            else:
                raise Exception(
                    f"Gene {gene.id} CDSs and exons under the mRNA do not match (different count)"
                )

        print(f"Gene {gene.id}: move {len(cdss)} CDSs to the mRNA")
        # No more issues? move the CDSs
        mrna.sub_features += cdss
        # And remove them from the gene
        gene.sub_features = gene_subf_clean
        gene.sub_features.append(mrna)

        return gene

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
                    raise Exception(f"Can't remove extra exons for {gene.id}, not all start with 'id-'")

        return gene

    def gene_to_exon(self, gene: SeqFeature) -> SeqFeature:
        """Returns an intermediary transcript for a gene with direct exon children."""

        transcript = SeqFeature(gene.location, type="mRNA")
        transcript.qualifiers["source"] = gene.qualifiers["source"]
        transcript.sub_features = []

        for exon in gene.sub_features:
            transcript.sub_features.append(exon)

        return transcript

    def add_funcann_feature(self, funcann: List, feature: SeqFeature, feat_type: str) -> None:
        """Append a feature object following the specifications.

        Args:
            funcann: The list to which the feature object are appended.
            feature: The SeqFeature to add to the list.
            feat_type: Feature type of the feature to store (e.g. gene, transcript, translation).

        """

        feature_object = {"object_type": feat_type, "id": feature.id}

        # Description?
        if "product" in feature.qualifiers:
            description = feature.qualifiers["product"][0]
            if self.check_product(description):
                feature_object["description"] = description

        if "Name" in feature.qualifiers and "description" not in feature_object:
            name = feature.qualifiers["Name"][0]

            # Exclude Name if it just a variant of the feature ID
            if feature.id not in name:
                feature_object["description"] = name

        # Synonyms?
        if "Name" in feature.qualifiers:
            feat_name = feature.qualifiers["Name"][0]
            if feat_name != feature.id:
                feature_object["synonyms"] = {"synonym": feat_name, "default": True}

        # is_pseudogene?
        if feature.type.startswith("pseudogen"):
            feature_object["is_pseudogene"] = True

        funcann.append(feature_object)

    def normalize_gene_id(self, gene: SeqFeature) -> str:
        """Remove any unnecessary prefixes around the gene ID.

        Generate a new stable id if it is not recognized as valid.

        Args:
            gene: Gene feature to normalize.

        Returns:
            A normalized gene id.

        """

        prefixes = ("gene-", "gene:")
        new_gene_id = self.remove_prefixes(gene.id, prefixes)

        # In case the gene id is not valid, use the GeneID
        if not self.valid_id(new_gene_id):
            print("Gene id is not valid: %s" % new_gene_id)
            qual = gene.qualifiers
            if "Dbxref" in qual:
                for xref in qual["Dbxref"]:
                    (db, value) = xref.split(":")
                    if db == "GeneID":
                        new_gene_id = db + "_" + value
                        print(f"Using GeneID {new_gene_id} for stable_id instead of {gene.id}")
                        return new_gene_id

            # Make a new stable_id
            if self.make_missing_stable_id:
                new_id = self.generate_stable_id()
                print(f"New id: {new_gene_id} -> {new_id}")
                return new_id
            else:
                raise Exception(f"Can't use invalid gene id for {gene}")

        return new_gene_id

    def generate_stable_id(self) -> str:
        """Returns a new unique gene stable_id with a prefix.

        The id is made up of a prefix and a number, which is auto incremented.
        Define the prefix with the param "stable_id_prefix",
        or use the genome organism_abbrev and prepend "TMP_" to it.

        """
        genome_data = InputSchema().genome_data
        if self.stable_id_prefix:
            prefix = self.stable_id_prefix
        else:
            dat = json.loads(genome_data.read())
            org = dat["BRC4"]["organism_abbrev"]
            prefix = "TMP_" + org + "_"
            self.stable_id_prefix = prefix

        if self.current_stable_id_number:
            number = self.current_stable_id_number
        else:
            number = 1

        number += 1
        new_id = "%s%d" % (prefix, number)
        self.current_stable_id_number = number

        return new_id

    def valid_id(self, name: str) -> bool:
        """Check that the format of a stable id is valid."""

        if not self.validate_gene_id:
            return True

        min_length = self.min_id_length

        # Trna (from tRNAscan)
        if re.search(r"^Trna", name):
            print("Stable id is a Trna from tRNA-scan: %s" % name)
            return False

        # Coordinates
        elif re.search(r"^.+:\d+..\d+", name):
            print("Stable id is a coordinate: %s" % name)
            return False

        # Special characters
        elif re.search(r"[ |]", name):
            print("Stable id contains special characters: %s" % name)
            return False

        # Min length
        elif len(name) <= min_length:
            print("Stable id is too short (<%d) %s" % (min_length, name))
            return False
        else:
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

        prefixes = ("cds-", "cds:")
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
                    if feat.id == "" or gene.id == feat.id:
                        feat.id = "%s_cds" % transcript.id
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


def main() -> None:
    """Main script entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    in_gff_path = mod.args["in_gff_path"]
    in_gff_path = Path(in_gff_path)
    out_gff_path = Path("gene_models.gff3")
    out_funcann_path = Path("functional_annotation.json")

    # Merge multiline gene features
    interim_gff_fh = tempfile.TemporaryFile(mode="w+")
    gff = process_gff3()
    gff.merge_genes_gff(in_gff_path, interim_gff_fh)
    interim_gff_fh.seek(0)
    # Load gff3 data and write a simpler version that follows our specifications
    gff.simpler_gff3(interim_gff_fh, out_gff_path, out_funcann_path)


if __name__ == "__main__":
    main()
