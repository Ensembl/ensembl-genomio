#!env python3

import eHive
from Bio import SeqIO
from os import path, rename
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
            "pseudogenic_transcript",
            "pseudogene",
            ]
    ignored_biotypes = [
            "biological_region",
            "supercontig",
            "chromosome",
            "scaffold",
            "contig",
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
            # Ensembl adds a prefix to features because different
            # types can use the same id in EG (gene and transcript)
            ensembl_mode : False
        }

    def run(self):
        gff3_path = self.param_required("gff_file")
        new_gff3_path = self.get_gff3(gff3_path)

        self.dataflow({ "specifications_gff_file" : new_gff3_path }, 2)

    def get_gff3(self, gff3_path):

        gff3_tmp_path = gff3_path + ".specifications.gff3"

        if gff3_path.endswith(".gz"):
            gff3_tmp_path += ".gz"
            with io.TextIOWrapper(gzip.open(gff3_path, "r")) as gff3_handle:
                with io.TextIOWrapper(gzip.open(gff3_tmp_path, "w")) as tmph:
                    self.parse_gff3(gff3_handle, tmph)
        else:
            with open(gff3_path, "r") as gff3_handle:
                with open(gff3_tmp_path, "w") as tmph:
                    self.parse_gff3(gff3_handle, tmph)


        return gff3_tmp_path

    def parse_gff3(self, gff3_handle, tmph):
        remove_prefix = True
        if self.param("ensembl_mode"):
            remove_prefix = False

        skipped_known_biotypes = {}
        skipped_unknown_biotypes = {}
        skipped_attributes = {}
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
                        if biotype in self.ignored_biotypes:
                            if biotype not in skipped_known_biotypes:
                                skipped_known_biotypes[biotype] = 0
                            skipped_known_biotypes[biotype] += 1
                        else:
                            if biotype not in skipped_unknown_biotypes:
                                skipped_unknown_biotypes[biotype] = 0
                            skipped_unknown_biotypes[biotype] += 1
                        continue

                # Filter attributes
                attrs = attribs.split(";")
                attributes = []
                for a in attrs:
                    (key, value) = a.split("=")
                    if (self.allowed_attributes and key not in self.allowed_attributes) or (biotype.endswith("gene") and key == 'Name'):
                        if not key in skipped_attributes:
                            skipped_attributes[key] = 0
                        skipped_attributes[key] += 1
                        continue

                    # Special production-specific modification
                    if remove_prefix and key in ("ID", "Parent") and ":" in value:

                        (btype, real_value) = value.split(":", 1)
                        if btype and real_value:
                            value = real_value


                    attributes.append(key + "=" + value)

                # Print back the line
                new_attribs = ";".join(attributes)
                new_line = "\t".join([chrom, source, biotype, start, end, col6, col7, col8, new_attribs]) + "\n"
                tmph.write(new_line)

        # Print stats
        if skipped_attributes:
            print("Skipped attributes:")
            for at in skipped_attributes:
                print("\t%s: %d" % (at, skipped_attributes[at]))

        if skipped_known_biotypes:
            print("Skipped known biotypes:")
            for bt in skipped_known_biotypes:
                print("\t%s: %d" % (bt, skipped_known_biotypes[bt]))
        if skipped_unknown_biotypes:
            print("Skipped unknown biotypes:")
            for bt in skipped_unknown_biotypes:
                print("\t%s: %d" % (bt, skipped_unknown_biotypes[bt]))
            raise Exception("Encountered an unknown biotype: %s. You need to say if it should be kept or ignored" % bt)

