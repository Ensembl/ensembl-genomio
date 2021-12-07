#!env python3

import re, sys, argparse
import os, json

from Bio import SeqIO
import tempfile

class FormattedFilesGenerator():
    """
    Contains a parser to load data from a file, and output a set of files that follow our schema for input into a core database
    """
    def __init__(self, output_dir):
        self.output_dir = output_dir
        try:
            os.mkdir(self.output_dir)
        except FileExistsError:
            pass
            
        
        self.dna_fasta = os.path.join(self.output_dir, 'dna.fasta')
        self.pep_fasta = os.path.join(self.output_dir, 'pep.fasta')
        self.genes_gff = os.path.join(self.output_dir, 'genes.gff')
        self.genome_json = os.path.join(self.output_dir, 'genome.json')
        self.seq_regions_json = os.path.join(self.output_dir, 'seq_regions.json')

        self.seqs = []
        
    def parse_genbank(self, gb_file):
        """
        Load a sequence from a Genbank file
        """
        
        with open(gb_file, "r") as gbh:
            for record in SeqIO.parse(gbh, "genbank"):
                # We don't want the record description (especially for the fasta file)
                record.description = ""
                self.seqs.append(record)
            
            self._write_dna_fasta()
            self._write_seq_regions_json()
    
    def _write_dna_fasta(self):
        with open(self.dna_fasta, "w") as fasta_fh:
            SeqIO.write(self.seqs, fasta_fh, "fasta")
    
    def _write_seq_regions_json(self):
        json_array = []
        for seq in self.seqs:
            seq_obj = {
                    "name" : seq.id,
                    "coord_system_level" : "chromosome",
                    "circular" : False,
                    "codon_table" : 1,
                    "length" : len(seq.seq)
                    }
            json_array.append(seq_obj)
        with open(self.seq_regions_json, "w") as seq_fh:
            seq_fh.write(json.dumps(json_array, indent=4))

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Prepare a genbank file for an additional sequence')
    
    parser.add_argument('--gb', type=str, required=True,
                help='Genbank file')
    parser.add_argument('--output_dir', type=str, required=True,
                help='Output dir')


    args = parser.parse_args()
    formatter = FormattedFilesGenerator(args.output_dir)
    formatter.parse_genbank(args.gb)

if __name__ == "__main__":
    main()
