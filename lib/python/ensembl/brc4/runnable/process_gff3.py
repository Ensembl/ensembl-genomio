#!env python3

import os, re, shutil
import eHive
import gzip
import csv, json
import tempfile

from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

class process_gff3(eHive.BaseRunnable):

    def param_defaults(self):
        return {
                "gene_types" : ("gene", "ncRNA_gene", "pseudogene"),
                "transcript_types" : (
                    "transcript",
                    "mRNA",
                    "pseudogenic_transcript",
                    "tRNA",
                    "pseudogenic_tRNA",
                    "rRNA",
                    "lnc_RNA",
                    "snoRNA",
                    "snRNA",
                    "ncRNA",
                    ),
                "ignored_types" : (
                    "region",
                    "gap",
                    "sequence_feature",
                    "repeat_region",
                    "mobile_genetic_element",
                    "microsatellite",
                    "cDNA_match"
                    ),
                "ncRNA_gene_types" : ("tRNA", "rRNA"),
                "skip_unrecognized" : False
                }
        

    def run(self):
        genome_data = self.param('genome_data')
        work_dir = self.param('work_dir')
        in_gff_path = self.param('in_gff3')

        # Create dedicated work dir
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)

        # Final files
        out_gff_path = os.path.join(work_dir, "gene_models.gff3")
        out_funcann_path = os.path.join(work_dir, "functional_annotation.json")

        # Merge multiline gene features
        interm_gff = tempfile.TemporaryFile(mode = "w+")
        self.merge_genes_gff(in_gff_path, interm_gff)
        interm_gff.seek(0)

        # Load gff3 data and write a simpler version that follows our specifications
        self.simpler_gff3(interm_gff, out_gff_path, out_funcann_path)

        # Output the gff3 file
        output = {
                "gff3": out_gff_path
                }
        self.dataflow(output, 2)
        
        # Output the functional annotation file
        output = {
                "metadata_type" : "functional_annotation",
                "metadata_json": out_funcann_path
                }
        self.dataflow(output, 3)

    def merge_genes_gff(self, in_gff_path, out_gff) -> None:
        """
        Merge genes in a gff that are split in multiple lines
        """
        tomerge = []
        
        with open(in_gff_path, "r") as gff3_in:
            for line in gff3_in:

                # Skip comments
                if line.startswith("#"):
                    out_gff.write(line)
                else:
                    line = line.rstrip()
                    fields = line.split("\t")
                    attr_fields = fields[8].split(";")
                    attrs = {}
                    for a in attr_fields:
                        (key, value) = a.split("=")
                        attrs[key] = value
                        
                    # Check this is a gene to merge; cache it then
                    if fields[2] == "gene" and "part" in attrs:
                        tomerge.append(fields)
                    
                    # If not, merge if needed, and print the line
                    else:
                        if tomerge:
                            for l in tomerge: print("\t".join(l))
                            new_line = self.merge_genes(tomerge)
                            print(new_line)
                            out_gff.write(new_line)
                            tomerge = []
                        out_gff.write(line + "\n")
                    
            # Print last merged gene if there is one
            if tomerge:
                new_line = self.merge_genes(tomerge)
                out_gff.write(new_line)


    def merge_genes(self, tomerge) -> str:
        
        min_start = -1
        max_end = -1
        for gene in tomerge:
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
                

    def simpler_gff3(self, gff3_in, out_gff_path, out_funcann_path) -> None:
        """
        Load a GFF3 from INSDC, and rewrite it in a simpler version,
        and also write a functional_annotation file
        """
        
        allowed_gene_types = self.param("gene_types")
        allowed_transcript_types = self.param("transcript_types")
        ignored_types = self.param("ignored_types")
        ncRNA_gene_types = self.param("ncRNA_gene_types")
        skip_unrecognized = self.param("skip_unrecognized")
        
        functional_annotation = []
        
        with open(out_gff_path, "w") as gff3_out:
            gff = GFF.parse(gff3_in)
            
            new_records = []
            fail_types = {}
            
            for record in gff:
                new_record = SeqRecord(record.seq, id=record.id)
                
                # GENES
                for gene in record.features:
                    
                    if gene.type in ignored_types:
                        continue
                    
                    if gene.type in ncRNA_gene_types:
                        # Transcript-level gene: add a gene parent
                        gene = self.ncrna_gene(gene)
                        
                    if gene.type in allowed_gene_types:
                        
                        # New gene ID 
                        gene.id = self.normalize_gene_id(gene.id)
                        
                        # Store gene functional annotation
                        self.add_funcann_feature(functional_annotation, gene, "gene")
                        
                        # replace qualifiers
                        old_qualifiers = gene.qualifiers
                        gene.qualifiers = {
                                "ID" : gene.id,
                                "source" : old_qualifiers["source"]
                                }
                        
                        # Gene with no subfeatures: need to create a transcript at least
                        if len(gene.sub_features) == 0:
                            print("Insert transcript for lone gene %s" % (gene.id))
                            transcript = self.transcript_for_gene(gene)
                            gene.sub_features = [transcript]
                        
                        # Transform gene - CDS to gene-transcript-exon-CDS
                        if gene.sub_features[0].type == "CDS":
                            print("Insert transcript-exon for %s (%d CDSs)" % (gene.id, len(gene.sub_features)))
                            transcript = self.gene_to_cds(gene)
                            gene.sub_features = [transcript]

                        # Transform gene - exon to gene-transcript-exon
                        if gene.sub_features[0].type == "exon":
                            print("Insert transcript for %s (%d exons)" % (gene.id, len(gene.sub_features)))
                            transcript = self.gene_to_exon(gene)
                            gene.sub_features = [transcript]

                        # TRANSCRIPTS
                        transcripts_to_delete = []
                        for count, transcript in enumerate(gene.sub_features):

                            if transcript.type not in allowed_transcript_types:
                                fail_types[transcript.type] = 1
                                message = "Unrecognized transcript type: %s for %s (%s)" % (transcript.type, transcript.id, gene.id)
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
                            old_qualifiers = transcript.qualifiers
                            transcript.qualifiers = {
                                    "ID" : transcript.id,
                                    "Parent" : gene.id,
                                    }
                            if "source" in old_qualifiers:
                                transcript.qualifiers["source"] = old_qualifiers["source"]

                            # EXONS AND CDS
                            cds_found = False
                            exons_to_delete = []
                            for tcount, feat in enumerate(transcript.sub_features):
                                
                                if feat.type == "exon":
                                    # Replace qualifiers
                                    old_qualifiers = feat.qualifiers
                                    feat.qualifiers = {
                                            "Parent" : transcript.id,
                                            }
                                    if "source" in old_qualifiers:
                                        feat.qualifiers["source"] = old_qualifiers["source"]
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
                                            "ID" : feat.id,
                                            "Parent" : transcript.id,
                                            "phase" : feat.qualifiers["phase"],
                                            "source" : feat.qualifiers["source"]
                                            }
                                else:
                                    fail_types[feat.type] = 1
                                    message = "Unrecognized exon type: %s" % exon.type
                                    print(message)
                                    if skip_unrecognized:
                                        exons_to_delete.append(tcount)
                                        continue
                            
                            if exons_to_delete:
                                for elt in sorted(exons_to_delete, reverse=True):
                                    transcript.sub_features.pop(elt)
                        
                        if transcripts_to_delete:
                            for elt in sorted(transcripts_to_delete, reverse=True):
                                gene.sub_features.pop(elt)

                        
                        # PSEUDOGENE CDS IDs
                        if gene.type == "pseudogene":
                            self.normalize_pseudogene_cds(gene)
                                
                    else:
                        fail_types[gene.type] = 1
                        message = "Unrecognized gene type: %s" % gene.type
                        print(message)
                        if skip_unrecognized:
                            del gene
                            continue

                    new_record.features.append(gene)
                new_records.append(new_record)
            
            if fail_types and not skip_unrecognized: raise Exception("Unrecognized types found (%s): fail" % (" ".join(fail_types.keys())))
            
            GFF.write(new_records, gff3_out)
        
        # Write functional annotation
        self.print_json(out_funcann_path, functional_annotation)
    
    def ncrna_gene(self, ncrna):
        """Create a gene for ncRNAs"""
        
        gene = SeqFeature(ncrna.location, type="ncRNA_gene")
        gene.qualifiers["source"] = ncrna.qualifiers["source"]
        gene.sub_features = [ncrna]
        gene.id = ncrna.id

        return gene
        

    def transcript_for_gene(self, gene):
        """Create a transcript for a lone gene"""
        
        transcript = SeqFeature(gene.location, type="mRNA")
        transcript.qualifiers["source"] = gene.qualifiers["source"]
        transcript.sub_features = []
        
        return transcript
        
    
    def gene_to_cds(self, gene):
        """Create a transcript - exon - cds chain"""
        
        transcript = SeqFeature(gene.location, type="mRNA")
        transcript.qualifiers["source"] = gene.qualifiers["source"]
        transcript.sub_features = []

        for cds in gene.sub_features:
            exon = SeqFeature(cds.location, type="exon")
            exon.qualifiers["source"] = gene.qualifiers["source"]
            transcript.sub_features.append(exon)
            transcript.sub_features.append(cds)
        
        return transcript
        
    def gene_to_exon(self, gene):
        """Create a transcript - exon chain"""
        
        transcript = SeqFeature(gene.location, type="mRNA")
        transcript.qualifiers["source"] = gene.qualifiers["source"]
        transcript.sub_features = []

        for exon in gene.sub_features:
            transcript.sub_features.append(exon)
        
        return transcript
        
    
    def print_json(self, path, data) -> None:
        """Dump an object to a json file"""
        
        with open(path, "w") as json_out:
            json_out.write(json.dumps(data, sort_keys=True, indent=4))
    
    def add_funcann_feature(self, funcann, feature, name):
        """Append a feature object following the specifications"""
        
        feature_object = {
                "object_type" : name,
                "id" : feature.id
                }
        
        # Description?
        if "product" in feature.qualifiers:
            description = feature.qualifiers["product"][0]
            if not re.search("^hypothetical protein$", description):
                feature_object["description"] = description
        
        # Synonyms?
        
        # is_pseudogene?
        if feature.type.startswith("pseudogen"):
            feature_object["is_pseudogene"] = True
        
        funcann.append(feature_object)
    
    def normalize_gene_id(self, gene_id) -> str:
        """Remove any unnecessary prefixes around the gene ID"""

        prefixes = ("gene-", "gene:")
        gene_id = self.remove_prefixes(gene_id, prefixes)
        return gene_id
    
    def normalize_transcript_id(self, gene_id, number) -> str:
        """Use a gene ID and a number to make a formatted transcript ID"""

        transcript_id = "%s_t%d" % (gene_id, number)
        return transcript_id

    def normalize_cds_id(self, cds_id) -> str:
        """
        Check the CDS ID is proper:
        - Remove any unnecessary prefixes around the CDS ID
        - Delete the ID if it is not proper
        """

        prefixes = ("cds-", "cds:")
        cds_id = self.remove_prefixes(cds_id, prefixes)

        # Special case: if the ID doesn't look like one, remove it
        if re.match("^...\|", cds_id):
            cds_id = ""
        
        return cds_id
    
    def normalize_pseudogene_cds(self, gene):
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
    
    def make_transcript_id(self, gene_id, transcript_number) -> str:
        """Create a transcript ID based on a gene and the number of the transcript"""
        
        # Simply add a numbered suffix to the gene_id
        transcript_id = "%s_t%d" (gene_id, transcript_number)

        return transcript_id

    def remove_prefixes(self, identifier, prefixes) -> str:
        """
        Remove prefixes from an identifier if they are found
        Return the unaltered identifier otherwise
        """
        for prefix in prefixes:
            if identifier.startswith(prefix):
                identifier = identifier[len(prefix):]
        return identifier

