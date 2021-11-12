#!env python3

import eHive
import gzip
import json
import io
import re
import sys
from functools import partial

from Bio import SeqIO
from os import path


class compare_report(eHive.BaseRunnable):

    def param_defaults(self):
        return {
        }

    def run(self):
        stats = self.param("stats")
        output_dir = self.param_required("output_dir")
        
        # Print report
        report = output_dir + "/report.log"
        print("Write report in %s" % report)
        
        fields = ("species", "accession", "seq_count_1", "seq_count_2", "num_diff_seq", "common", "only1", "only2", "max_only1", "max_only2","other_locations","summary","organellar_summary")
        
        with open(report, "w") as out_fh:
            out_fh.write("#" + "\t".join(fields) + "\n")
            
            for species in sorted(stats.keys()):
                stat = stats[species]
                stat["species"] = species
                
                line = []
                for f in fields:
                    if f in stat:
                        line.append(str(stat[f]))
                out_fh.write("\t".join(line) + "\n")
                
