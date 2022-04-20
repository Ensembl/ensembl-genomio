#!env python3

import re, sys, argparse
import os, json
import eHive

from Bio import SeqIO
from Bio import GenBank

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF

import tempfile
from collections import Counter

class extract_from_gb(eHive.BaseRunnable):

    def param_defaults(self):
        return {
                }

    def run(self):
        accession = self.param_required('gb_accession')
        work_dir = self.param('work_dir')
        prefix = self.param('ids_prefix')
        prod_name = self.param('production_name')
        gb_file = self.param('gb_file')

        formatter = FormattedFilesGenerator(work_dir)
        formatter.set_prefix(prefix)
        formatter.set_production_name(prod_name)
        formatter.parse_genbank(gb_file)

        # Output the gff3 file
        output = {
                "gff3": formatter.genes_gff,
                "fasta_dna": formatter.fasta_dna,
                "fasta_pep": formatter.fasta_pep,
                "seq_region": formatter.seq_regions,
                "genome_data": formatter.genome,
                }
        self.dataflow(output, 2)

class FormattedFilesGenerator():
    """
    Contains a parser to load data from a file, and output a set of files that follow our schema for input into a core database
    """
    
    locations = {
            "mitochondrion" : "mitochondrial_chromosome",
            "apicoplast" : "apicoplast_chromosome",
            "chloroplast" : "chloroplast_chromosome",
            "chromoplast" : "chromoplast_chromosome",
            "cyanelle" : "cyanelle_chromosome",
            "leucoplast" : "leucoplast_chromosome",
            }

    allowed_feat_types = [
            'gene',
            'transcript',
            'tRNA',
            'rRNA',
            'CDS',
            ]

    
    def __init__(self, output_dir, prefix=""):
        self.output_dir = output_dir
        try:
            os.makedirs(self.output_dir)
        except FileExistsError:
            pass
            
        self.genome = os.path.join(self.output_dir, 'genome.json')
        self.seq_regions = os.path.join(self.output_dir, 'seq_regions.json')
        self.fasta_dna = os.path.join(self.output_dir, 'dna.fasta')
        self.fasta_pep = os.path.join(self.output_dir, 'pep.fasta')
        self.genes_gff = os.path.join(self.output_dir, 'genes.gff')
        self.prefix = prefix
        self.seqs = []
    
    def set_prefix(self, prefix):
        """
        Define a prefix to add to the feature IDs
        """
        if prefix:
            self.prefix = prefix
    
    def set_production_name(self, prod_name):
        """
        Define a production_name for the genome
        """
        if prod_name:
            self.prod_name = prod_name
        
    def parse_genbank(self, gb_file):
        """
        Load a sequence from a Genbank file
        """
        
        organella = self._get_organella(gb_file)
        
        with open(gb_file, "r") as gbh:
            for record in SeqIO.parse(gbh, "genbank"):
                # We don't want the record description (especially for the fasta file)
                record.description = ""
                record.organelle = None
                if record.id in organella:
                    record.organelle = organella[record.id]
                self.seqs.append(record)
            
            self._write_genome_json()
            self._write_genes_gff()
            self._write_seq_regions_json()
            self._write_fasta_dna()
    
    def _get_organella(self, gb_file):
        """
        Retrive the organelle from the genbank file, using the specific GenBank object,
        because SeqIO does not support this field
        """
        organella = {}
        with open(gb_file, "r") as gbh:
            for record in GenBank.parse(gbh):
                accession = record.version
                for q in record.features[0].qualifiers:
                    if q.key == "/organelle=":
                        organelle = q.value.replace('"', '')
                        organella[record.version] = organelle
        return organella
    
    def _write_fasta_dna(self):
        with open(self.fasta_dna, "w") as fasta_fh:
            SeqIO.write(self.seqs, fasta_fh, "fasta")
    
    def _write_genes_gff(self):
        peptides = []
        
        with open(self.genes_gff, "w") as gff_fh:
            recs = []
            all_ids = []
            
            for seq in self.seqs:
                feats = {}
                
                for feat in seq.features:
                    if feat.type not in self.allowed_feat_types: continue
                    gff_qualifiers = feat.qualifiers
                    gff_feat = SeqFeature(
                            location=feat.location,
                            type=feat.type,
                            strand=feat.location.strand,
                            qualifiers=gff_qualifiers
                            )
                    
                    if "gene" in gff_qualifiers:
                        gene_name = gff_qualifiers["gene"][0]
                        gene_id = self.prefix + gene_name
                        
                        if feat.type == "gene":
                            gff_feat.qualifiers["ID"] = gene_id
                            gff_feat.qualifiers["Name"] = gene_name
                            del gff_feat.qualifiers["gene"]
                            feats[str(gene_id)] = gff_feat
                            all_ids.append(str(gene_id))
                        
                        if feat.type == "CDS":
                            cds_id = gene_id + "_p1"
                            tr_id = gene_id + "_t1"
                            gff_feat.qualifiers["ID"] = cds_id
                            gff_feat.qualifiers["Parent"] = tr_id
                            del gff_feat.qualifiers["gene"]
                            
                            # Add fasta to pep fasta file
                            peptides.append(
                                SeqRecord(
                                    Seq(feat.qualifiers["translation"][0]),
                                    id=cds_id
                                    )
                                )

                            # Also create a parent transcript for this translation
                            tr_qualifiers = {
                                    "ID" : tr_id,
                                    "Name" : gene_name,
                                    "Parent": gene_id
                                    }
                            gff_tr = SeqFeature(
                                    location=feat.location,
                                    type="mRNA",
                                    strand=feat.location.strand,
                                    qualifiers=tr_qualifiers
                                    )
                            feats[str(tr_id)] = gff_tr
                            feats[str(cds_id)] = gff_feat
                            all_ids.append(str(tr_id))
                            all_ids.append(str(cds_id))
                        
                    elif feat.type in ("tRNA", "rRNA"):
                            feat_name = gff_qualifiers["product"][0]
                            gene_id = self.prefix + feat_name

                            parts = gene_id.split(" ")
                            if len(parts) > 2:
                                print(f"Shortening gene_id to {parts[0]}")
                                gene_id = parts[0]
                            gene_id = self._uniquify_id(gene_id, all_ids)

                            feat_id = gene_id + "_t1"
                            gff_feat.qualifiers["ID"] = feat_id
                            gff_feat.qualifiers["Name"] = feat_name
                            gff_feat.qualifiers["Parent"] = gene_id
                            
                            # Also create a parent gene for this transcript
                            gene_qualifiers = {
                                    "ID" : gene_id,
                                    "Name" : feat_name,
                                    }
                            gff_gene = SeqFeature(
                                    location=feat.location,
                                    type="gene",
                                    strand=feat.location.strand,
                                    qualifiers=gene_qualifiers
                                    )
                            feats[str(gene_id)] = gff_gene
                            feats[str(feat_id)] = gff_feat
                            all_ids.append(str(gene_id))
                            all_ids.append(str(feat_id))
                            
                rec = SeqRecord(seq.seq, seq.id)
                rec.features = feats.values()
                recs.append(rec)
    
            GFF.write(recs, gff_fh)

            with open(self.fasta_pep, "w") as fasta_fh:
                SeqIO.write(peptides, fasta_fh, "fasta")
            
            # Warn if some IDs are not unique
            count = dict(Counter(all_ids))
            num_duplicates = 0
            for key in count:
                if count[key] > 1:
                    num_duplicates += 1
                    print(f"ID {key} is duplicated {count[key]} times")
            if num_duplicates > 0:
                raise Exception(f"Some {num_duplicates} IDs are duplicated")

    def _uniquify_id(self, gene_id, all_ids):
        """Ensure the gene id used is unique,
        and append a number otherwise, starting at 2
        """
        
        new_id = gene_id
        num = 1
        while new_id in all_ids:
            print(f"{new_id} exists, update")
            num += 1
            new_id = f"{gene_id}_{num}"
        print(f"Using {new_id}")
        
        return new_id
        
        
    def _write_seq_regions_json(self):
        json_array = []
        
        for seq in self.seqs:
            codon_table = self._get_codon_table(seq)
            if not codon_table:
                print(f"Warning: No codon table found. Make sure to change the codon table number in {self.seq_regions_json} manually if it is not the standard codon table")

                codon_table = 1
            else:
                codon_table = int(codon_table)
            seq_obj = {
                    "name" : seq.id,
                    "coord_system_level" : "chromosome",
                    "circular" : (seq.annotations["topology"] == "circular"),
                    "codon_table" : codon_table,
                    "length" : len(seq.seq),
                    }
            if seq.organelle:
                seq_obj["location"] = self._prepare_location(seq.organelle)
                if not codon_table:
                    print(f"Warning: '{seq.organelle}' is an organelle: make sure to change the codon table number in {self.seq_regions_json} manually if it is not the standard codon table")

            # Additional attributes for Ensembl
            seq_obj["added_sequence"] = {
                    "accession" : seq.id,
                    "assembly_provider" : {
                        "name" : "GenBank",
                        "url" : "https://www.ncbi.nlm.nih.gov/genbank",
                        }
                    }
            if not seq_obj["added_sequence"]["assembly_provider"]["name"]:
                print(f"Warning: please add the relevant provider name for the assembly in {self.seq_regions}")
            if not seq_obj["added_sequence"]["assembly_provider"]["url"]:
                print(f"Warning: please add the relevant provider url for the assembly in {self.seq_regions}")

            # Additional attributes for gene set, if any
            # TODO

            json_array.append(seq_obj)
        with open(self.seq_regions, "w") as seq_fh:
            seq_fh.write(json.dumps(json_array, indent=4))

    def _get_codon_table(self, seq):
        """
        Look at the CDS features to see if they have a codon table
        """
        for feat in seq.features:
            if feat.type == "CDS":
                quals = feat.qualifiers
                if "transl_table" in quals:
                    return quals["transl_table"][0]
                else:
                    return
        return

    def _prepare_location(self, organelle):
        """
        Given an organelle name, returns the SO term corresponding to its location
        """
        if organelle in self.locations:
            return self.locations[organelle]
        else:
            raise Exception(f"Unkown organelle: {location}")

    def _write_genome_json(self):
        """
        Write a draft for the genome json file
        Only the production_name is needed, but the rest of the fields need to be given
        for the validation of the json file
        """
        
        prod_name = self.prod_name if self.prod_name else ""
        
        genome_data = {
                "species" : {
                    "production_name" : prod_name,
                    "taxonomy_id" : 0,
                    },
                "assembly" : {
                    "accession" : "GCA_000000000",
                    "version" : 1
                    },
                "added_seq" : {}
                }

        if not genome_data["species"]["production_name"]:
            print(f"Warning: please add the relevant production_name for this genome in {self.genome}")
            
        ids = [ seq.id for seq in self.seqs]
        genome_data["added_seq"]["region_name"] = ids
        
        with open(self.genome, "w") as genome_fh:
            genome_fh.write(json.dumps(genome_data, indent=4))
