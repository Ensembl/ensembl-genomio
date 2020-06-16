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
                "metadata_type" : "json",
                "metadata_json": out_funcann_path
                }
        self.dataflow(output, 3)

    def simpler_gff3(self, in_gff_path, out_gff_path, out_funcann_path):
        """
        Load a GFF3 from INSDC, and rewrite it in a simpler version,
        and also write a functional_annotation file
        """
        
        with open(out_gff_path, "w") as gff3_out:
            with open(in_gff_path, "r") as gff3_in:
                gff = GFF.parse(gff3_in)
                
                new_records = []
                for record in gff:
                    new_record = SeqRecord(record.seq, id=record.id)
                    
                    # GENES
                    for gene in record.features:
                        if gene.type in ("region"):
                            continue
                        if gene.type in ("gene", "ncRNA_gene", "pseudogene"):
                            gene.id = self.normalize_gene_id(gene.id)
                            gene.qualifiers = {
                                    "ID" : gene.id
                                    }

                            # TRANSCRIPTS
                            for count, transcript in enumerate(gene.sub_features):
                                if transcript.type in ("transcript", "mRNA", "pseudogenic_transcript", "tRNA", "pseudogenic_tRNA"):
                                    transcript_number = count + 1
                                    transcript.id = self.normalize_transcript_id(gene.id, transcript_number)
                                    transcript.qualifiers = {
                                            "ID" : transcript.id,
                                            "Parent" : gene.id,
                                            }
                                else:
                                    raise Exception("Unrecognized transcript type: %s" % transcript.type)

                                # EXONS AND CDS
                                for count, feat in enumerate(transcript.sub_features):
                                    if feat.type in ("exon"):
                                        feat.qualifiers = {
                                                "Parent" : transcript.id,
                                                }
                                    elif feat.type in ("CDS"):
                                        feat.id = self.normalize_cds_id(feat.id)
                                        feat.qualifiers = {
                                                "ID" : feat.id,
                                                "Parent" : transcript.id,
                                                }
                                    else:
                                        raise Exception("Unrecognized exon type: %s" % exon.type)
                                    
                        else:
                            raise Exception("Unrecognized gene type: %s" % gene.type)
                    
                        new_record.features.append(gene)
                    new_records.append(new_record)
                
                GFF.write(new_records, gff3_out)
                    
    
    def normalize_gene_id(self, gene_id) -> str:
        prefixes = ("gene-", "gene:")
        for prefix in prefixes:
            gene_id = self.remove_prefix(gene_id, prefix)
        return gene_id
    
    def normalize_transcript_id(self, gene_id, number) -> str:
        transcript_id = "%s_t%d" % (gene_id, number)
        return transcript_id

    def normalize_cds_id(self, cds_id) -> str:
        prefixes = ("cds-", "cds:")
        for prefix in prefixes:
            cds_id = self.remove_prefix(cds_id, prefix)
        return cds_id
    
    def make_transcript_id(self, gene_id, transcript_number) -> str:
        """
        Create a transcript ID based on a gene and the number of the transcript
        """
        # Simply add a numbered suffix to the gene_id
        transcript_id = "%s_t%d" (gene_id, transcript_number)

        return transcript_id

    def remove_prefix(self, identifier, prefix) -> str:
        """
        Remove a prefix from an identifier if it is found
        Return the unaltered identifier otherwise
        """
        if identifier.startswith(prefix):
            identifier = identifier[len(prefix):]
        return identifier

