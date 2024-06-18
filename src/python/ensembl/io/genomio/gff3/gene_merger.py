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
"""Merge split genes in a GFF file."""

__all__ = [
    "GFFGeneMerger",
]

from importlib.resources import as_file, files
import logging
from os import PathLike
from pathlib import Path
import re
from typing import List

import ensembl.io.genomio.data.gff3
from ensembl.io.genomio.utils.json_utils import get_json


class GFFGeneMerger:
    """Specialized class to merge split genes in a GFF3 file, prior to further parsing."""

    def __init__(self) -> None:
        source = files(ensembl.io.genomio.data.gff3).joinpath("biotypes.json")
        with as_file(source) as biotypes_json:
            self._biotypes = get_json(biotypes_json)

    def merge(self, in_gff_path: PathLike, out_gff_path: PathLike) -> List[str]:
        """
        Merge genes in a gff that are split in multiple lines.

        Args:
            in_gff_path: Input GFF3 that may have split merge.
            out_gff_path: Output GFF3 with those genes merged.

        Returns:
            List of all merged genes, each represented as a string of the GFF3 lines of all their parts.
        """
        to_merge = []
        merged: List[str] = []

        with Path(in_gff_path).open("r") as in_gff_fh, Path(out_gff_path).open("w") as out_gff_fh:
            for line in in_gff_fh:
                # Skip comments
                if line.startswith("#"):
                    if line.startswith("##FASTA"):
                        logging.warning("This GFF3 file contains FASTA sequences")
                        break
                    out_gff_fh.write(line)
                else:
                    # Parse one line
                    line = line.rstrip()
                    fields = line.split("\t")
                    attr_fields = fields[8].split(";")
                    attrs = {}
                    for a in attr_fields:
                        (key, value) = a.split("=")
                        attrs[key] = value

                    # Check this is a gene to merge; cache it then
                    if fields[2] in self._biotypes["gene"]["supported"] and (
                        "part" in attrs or "is_ordered" in attrs
                    ):
                        to_merge.append(fields)

                    # If not, merge previous gene if needed, and print the line
                    else:
                        if to_merge:
                            merged_str = []
                            for line_to_merge in to_merge:
                                merged_str.append("\t".join(line_to_merge))
                            merged.append("\n".join(merged_str) + "\n")

                            new_line = self._merge_genes(to_merge)
                            out_gff_fh.write(new_line)
                            to_merge = []
                        out_gff_fh.write(line + "\n")

            # Print last merged gene if there is one
            if to_merge:
                merged_str = []
                for line_to_merge in to_merge:
                    merged_str.append("\t".join(line_to_merge))
                merged.append("\n".join(merged_str) + "\n")

                new_line = self._merge_genes(to_merge)
                out_gff_fh.write(new_line)

        logging.debug(f"Merged lines: {len(merged)}")
        return merged

    def _merge_genes(self, to_merge: List) -> str:
        """Returns a single gene gff3 line merged from separate parts.

        Args:
            to_merge: List of gff3 lines with gene parts.

        """
        min_start = -1
        max_end = -1
        for gene in to_merge:
            start = int(gene[3])
            end = int(gene[4])

            if start < min_start or min_start < 0:
                min_start = start
            if end > max_end or max_end < 0:
                max_end = end

        # Take the first line as template and replace things
        new_gene = to_merge[0]
        new_gene[3] = str(min_start)
        new_gene[4] = str(max_end)

        attrs = new_gene[8]
        attrs = attrs.replace(";is_ordered=true", "")
        attrs = re.sub(r";part=\d+/\d+", "", attrs)
        new_gene[8] = attrs

        return "\t".join(new_gene) + "\n"
