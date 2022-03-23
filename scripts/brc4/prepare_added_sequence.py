#!env python3

import re, sys, argparse
import os, json
import hashlib

from Bio import SeqIO
from Bio import GenBank

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF

import tempfile

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
            os.mkdir(self.output_dir)
        except FileExistsError:
            pass
            
        
        self.manifest = os.path.join(self.output_dir, 'manifest.json')
        self.genome_json = os.path.join(self.output_dir, 'genome.json')
        self.seq_regions_json = os.path.join(self.output_dir, 'seq_regions.json')
        self.functional_annotation_json = os.path.join(self.output_dir, 'functional_annotation.json')
        self.dna_fasta = os.path.join(self.output_dir, 'dna.fasta')
        self.pep_fasta = os.path.join(self.output_dir, 'pep.fasta')
        self.genes_gff = os.path.join(self.output_dir, 'genes.gff')
        self.prefix = prefix
        self.seqs = []
    
    def set_prefix(self, prefix):
        """
        Define a prefix to add to the feature IDs
        """
        if prefix:
            self.prefix = prefix
        
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
            self._write_dna_fasta()
            self._write_manifest()
    
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
    
    def _write_dna_fasta(self):
        with open(self.dna_fasta, "w") as fasta_fh:
            SeqIO.write(self.seqs, fasta_fh, "fasta")
    
    def _write_genes_gff(self):
        funcs = []
        peptides = []
        
        with open(self.genes_gff, "w") as gff_fh:
            recs = []
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
                            
                            gene_ann = {
                                    "id" : gene_id,
                                    "type" : "gene",
                                    "synonyms" : [
                                        { "default" : True, "id" : gene_name }
                                        ]
                                    }
                            funcs.append(gene_ann)
                        
                        if feat.type == "CDS":
                            cds_id = gene_id + "_p1"
                            tr_id = gene_id + "_t1"
                            gff_feat.qualifiers["ID"] = cds_id
                            gff_feat.qualifiers["Parent"] = tr_id
                            gff_feat.qualifiers["Name"] = gene_name
                            del gff_feat.qualifiers["gene"]

                            trl_ann = {
                                    "id" : cds_id,
                                    "type" : "translation",
                                    }
                            funcs.append(trl_ann)
                            
                            # Add fasta to pep fasta file
                            peptides.append({ "id" : cds_id, "sequence": feat.qualifiers["translation"] })

                            # Also create a parent transcript for this translation
                            tr_ann = {
                                    "id" : tr_id,
                                    "type" : "transcript",
                                    "synonyms" : [
                                        { "default" : True, "id" : gene_name }
                                        ]
                                    }
                            funcs.append(tr_ann)

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
                        
                    elif feat.type in ("tRNA", "rRNA"):
                            feat_name = gff_qualifiers["product"][0]
                            gene_id = self.prefix + feat_name

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

                            gene_rna_ann = {
                                    "id" : gene_id,
                                    "type" : "gene",
                                    "description" : feat_name,
                                    "synonyms" : [
                                        { "default" : True, "id" : feat_name }
                                        ]
                                    }
                            funcs.append(gene_rna_ann)
                            rna_ann = {
                                    "id" : feat_id,
                                    "type" : "transcript",
                                    "description" : feat_name,
                                    }
                            funcs.append(rna_ann)
                            
                            
                
                rec = SeqRecord(Seq("X"), "Sequence")
                rec.features = feats.values()
                recs.append(rec)

            GFF.write(recs, gff_fh)
            
            with open(self.functional_annotation_json, "w") as func_json:
                func_json.write(json.dumps(funcs, indent=4))
            # TODO: fasta pep

            with open(self.fasta_pep, "w") as fpep:
                for seq in peptides:
                    fpep.write(">" + seq["id"] + "\n" + seq["seq"])
    
    def _write_manifest(self):
        """
        Write a manifest following the defined schema
        """
        manifest = {}
        self._add_to_manifest(manifest, 'gff3', self.genes_gff)
        self._add_to_manifest(manifest, 'fasta_dna', self.dna_fasta)
        self._add_to_manifest(manifest, 'fasta_pep', self.pep_fasta)
        self._add_to_manifest(manifest, 'seq_region', self.seq_regions_json)
        self._add_to_manifest(manifest, 'genome', self.genome_json)
        self._add_to_manifest(manifest, 'functional_annotation', self.functional_annotation_json)
        with open(self.manifest, "w") as manifest_fh:
            manifest_fh.write(json.dumps(manifest, indent=4))
    
    def _add_to_manifest(self, manifest, key, path):
        """
        Add a file to the manifest, if it exists
        """
        if os.path.exists(path):
            with open(path, "rb") as file_fh:
                content = file_fh.read()
                md5 = hashlib.md5(content).hexdigest()
                manifest[key] = {
                        "file": os.path.basename(path),
                        "md5sum": md5
                        }
        
    def _write_seq_regions_json(self):
        json_array = []
        
        for seq in self.seqs:
            codon_table = self._get_codon_table(seq)
            if not codon_table:
                print(f"Warning: No codon table found. Make sure to change the codon table number in {self.seq_regions_json} manually if it is not the standard codon table")

                codon_table = 1
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
                        "name" : "",
                        "url" : "",
                        }
                    }
            if not seq_obj["added_sequence"]["assembly_provider"]["name"]:
                print(f"Warning: please add the relevant provider name for the assembly in {self.seq_regions_json}")
            if not seq_obj["added_sequence"]["assembly_provider"]["url"]:
                print(f"Warning: please add the relevant provider url for the assembly in {self.seq_regions_json}")

            # Additional attributes for gene set, if any
            # TODO

            json_array.append(seq_obj)
        with open(self.seq_regions_json, "w") as seq_fh:
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
        genome_data = {
                "species" : {
                    "production_name" : "",
                    "taxonomy_id" : 0,
                    },
                "assembly" : {
                    "accession" : "GCA_000000000",
                    "version" : 1
                    },
                "added_seq" : {}
                }

        if not genome_data["species"]["production_name"]:
            print(f"Warning: please add the relevant production_name for this genome in {self.genome_json}")
            
        ids = [ seq.id for seq in self.seqs]
        genome_data["added_seq"]["region_name"] = ids
        
        with open(self.genome_json, "w") as genome_fh:
            genome_fh.write(json.dumps(genome_data, indent=4))

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Prepare a genbank file for an additional sequence')
    
    parser.add_argument('--gb', type=str, required=True,
                help='Genbank file')
    parser.add_argument('--output_dir', type=str, required=True,
                help='Output dir')
    parser.add_argument('--prefix', type=str,
                help='IDs prefix')

    args = parser.parse_args()
    formatter = FormattedFilesGenerator(args.output_dir)
    formatter.set_prefix(args.prefix)
    formatter.parse_genbank(args.gb)

if __name__ == "__main__":
    main()
