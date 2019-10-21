#!env python3

import eHive
from Bio import SeqIO
from os import path
import gzip
import io
from math import floor

class GFF3Specifier(eHive.BaseRunnable):

    allowed_biotypes = [
            "gene",
            "mRNA",
            "exon",
            "CDS",
            "three_prime_UTR",
            "five_prime_UTR",
            ]
    allowed_attributes = [
            "ID",
            "Parent",
            "Derives_from",
            "Name",
            "biotype",
            ]

    def param_default(self):
        return {
        }

    def run(self):
        gff3_path = self.param_required("gff3_file")

        gff = self.get_gff3(gff3_path)

    def get_gff3(self, gff3_path):

        seqs = {};
        genes = {};
        peps = {};

        gff3_tmp_path = gff3_path + ".tmp.gz"

        skipped_biotypes = {}
        skipped_attributes = {}

        with io.TextIOWrapper(gzip.open(gff3_path, "r")) as gff3_handle:
            with io.TextIOWrapper(gzip.open(gff3_tmp_path, "w")) as tmph:

                for line in gff3_handle:
                    # Skip separator lines
                    if line[0:3] == "###":
                        continue
                    # Keep all comment lines intact
                    elif line[0] == "#":
                        tmph.write(line)
                    else:
                        line = line.rstrip()
                        (chrom, source, biotype, start, end, col6, col7, col8, attribs) = line.split("\t")

                        # Filter biotype
                        if self.allowed_biotypes and biotype not in self.allowed_biotypes:
                            if not (biotype.endswith("RNA") or biotype.endswith("RNA_gene")):
                                if biotype not in skipped_biotypes:
                                    skipped_biotypes[biotype] = 0
                                skipped_biotypes[biotype] += 1
                                continue

                        # Filter attributes
                        attrs = attribs.split(";")
                        attributes = []
                        for a in attrs:
                            (key, value) = a.split("=")
                            if self.allowed_attributes and key not in self.allowed_attributes:
                                if not key in skipped_attributes:
                                    skipped_attributes[key] = 0
                                skipped_attributes[key] += 1
                                continue

                            # Special production-specific modification
                            if key == "ID":
                                (btype, real_value) = value.split(":")
                                value = real_value

                            attributes.append(key + ":" + value)

                        # Print back the line
                        new_attribs = ";".join(attributes)
                        new_line = "\t".join([chrom, source, biotype, start, end, col6, col7, col8, new_attribs]) + "\n"
                        tmph.write(new_line)

                        #raise Exception("die")

        if skipped_biotypes:
            print("Skipped biotypes:")
            for bt in skipped_biotypes:
                print("\t%s: %d" % (bt, skipped_biotypes[bt]))
        if skipped_attributes:
            print("Skipped attributes:")
            for at in skipped_attributes:
                print("\t%s: %d" % (at, skipped_attributes[at]))


        return gff3_tmp_path
