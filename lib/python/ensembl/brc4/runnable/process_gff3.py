#!env python3

import os, re, shutil
import eHive
import gzip
import csv, json

from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

class process_gff3(eHive.BaseRunnable):

    def param_defaults(self):
        return {
                "gene_types" : ("gene", "ncRNA_gene", "pseudogene"),
                "transcript_types" : ("transcript", "mRNA", "pseudogenic_transcript", "tRNA", "pseudogenic_tRNA", "rRNA"),
                "ignored_types" : ("region", "gap"),
                "ncRNA_gene_types" : ("tRNA", "rRNA"),
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

        # Load gff3 data and write a simpler version that follows our specifications
        self.simpler_gff3(in_gff_path, out_gff_path, out_funcann_path)

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

    def simpler_gff3(self, in_gff_path, out_gff_path, out_funcann_path) -> None:
        """
        Load a GFF3 from INSDC, and rewrite it in a simpler version,
        and also write a functional_annotation file
        """
        
        allowed_gene_types = self.param("gene_types")
        allowed_transcript_types = self.param("transcript_types")
        ignored_types = self.param("ignored_types")
        ncRNA_gene_types = self.param("ncRNA_gene_types")
        
        functional_annotation = []
        
        with open(out_gff_path, "w") as gff3_out:
            with open(in_gff_path, "r") as gff3_in:
                gff = GFF.parse(gff3_in)
                
                new_records = []
                for record in gff:
                    new_record = SeqRecord(record.seq, id=record.id)
                    
                    # GENES
                    for gene in record.features:
                        print(gene.type)
                        
                        if gene.type in ignored_types:
                            continue
                        
                        if gene.type in ncRNA_gene_types:
                            # Transcript-level gene: add a gene parent
                            gene = self.ncrna_gene(gene)
                            print(gene)
                            
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
                            
                            # Transform gene - CDS to gene-transcript-exon-CDS
                            if gene.sub_features[0].type == "CDS":
                                print("Insert transcript-exon for %s (%d CDSs)" % (gene.id, len(gene.sub_features)))
                                transcript = self.gene_to_cds(gene)
                                gene.sub_features = [transcript]
                            
                            # PSEUDOGENE CDS IDs
                            if gene.type == "pseudogene":
                                self.normalize_pseudogene_cds(gene)

                            # TRANSCRIPTS
                            for count, transcript in enumerate(gene.sub_features):

                                if transcript.type not in allowed_transcript_types:
                                    raise Exception("Unrecognized transcript type: %s for %s" % (transcript.type, transcript.id))

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
                                for count, feat in enumerate(transcript.sub_features):
                                    
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
                                        raise Exception("Unrecognized exon type: %s" % exon.type)
                                    
                        else:
                            print(gene)
                            raise Exception("Unrecognized gene type: %s" % gene.type)
                    
                        new_record.features.append(gene)
                    new_records.append(new_record)
                
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
        """Remove any unnecessary prefixes around the CDS ID"""

        prefixes = ("cds-", "cds:")
        cds_id = self.remove_prefixes(cds_id, prefixes)
        return cds_id
    
    def normalize_pseudogene_cds(self, gene):
        """Ensure CDS from a pseudogene have a different ID"""

        for transcript in gene.sub_features:
            for feat in transcript.sub_features:
                if feat.type == "CDS":
                    feat.id = self.normalize_cds_id(feat.id)
                    if gene.id == feat.id:
                        feat.id = "%s_cds" % gene.id
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

