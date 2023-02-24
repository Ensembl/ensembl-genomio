#!env python3
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


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
        return {}

    def run(self):
        stats = self.param("stats")
        output_dir = self.param_required("output_dir")

        # Print report
        report = output_dir + "/report.log"
        print("Write report in %s" % report)

        fields = (
            "species",
            "accession",
            "seq_count_1",
            "seq_count_2",
            "num_diff_seq",
            "common",
            "only1",
            "only2",
            "max_only1",
            "max_only2",
            "other_locations",
            "summary",
            "organellar_summary",
            "Assembly_level_1",
            "Assembly_level_2",
        )

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
