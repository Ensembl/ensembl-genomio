#!env python3

#from BCBio import GFF

import eHive
import json
from Bio import SeqIO
from os import path

class Integrity(eHive.BaseRunnable):

    def param_default(self):
        return {
        }

    def run(self):
        manifest_path = self.param_required('manifest')
        #print(manifest_path)

        errors = []

        with open(manifest_path) as manifest_file:
            manifest = json.load(manifest_file)
            
            # Use dir name from the manifest
            for f in manifest:
                manifest[f] = path.join(path.dirname(manifest_path), manifest[f])
            
            # Get content
            dna = {}
            pep = {}
            seq_regions = {}

            if "fasta_dna" in manifest:
                dna = self.get_fasta_lengths(manifest["fasta_dna"])
            if "fasta_pep" in manifest:
                pep = self.get_fasta_lengths(manifest["fasta_pep"])
            if "metadata_seq_region" in manifest:
                seq_regions = self.get_json(manifest["metadata_seq_region"])

            # Check fasta dna and seq_region integrity
            if dna and seq_regions:
                errors += self.check_dna_seq_regions(dna, seq_regions)

        if errors:
            raise Exception("Integrity test failed: " + str(errors))

    def get_fasta_lengths(self, fasta_path):
        data = {}
        for rec in SeqIO.parse(fasta_path, "fasta"):
            data[rec.id] = len(rec)
        return data

    def get_json(self, json_path):
        with open(json_path) as json_file:
            return json.load(json_file)

    def check_dna_seq_regions(self, dna, seqr):

        seqr_lengths = {}
        for seq in seqr:
            seqr_lengths[seq["name"]] = seq["length"]

        common = []
        diff = []
        only_dna = []
        only_seqr = []
        for s in dna:
            if s in seqr_lengths:
                dlength = dna[s]
                slength = seqr_lengths[s]

                if dlength == slength:
                    common.append(s)
                else:
                    diff.append(s)

            else:
                only_dna.append(s)

        for s in seqr_lengths:
            if s not in common and s not in diff:
                only_seqr.append(s)
        
        errors = []
        if common:
            print("%d seq_regions in common between fasta dna and seq_region" % len(common))
        if diff:
            errors.append("%d seq_regions are of different length between fasta dna and seq_region" % len(diff))
        if only_dna:
            errors.append("WARNING: %d seq_regions only found in fasta dna" % len(only_dna))
#            print(only_dna)
        if only_seqr:
            errors.append("%d seq_regions only found in seq_region" % len(only_seqr))

        return errors

