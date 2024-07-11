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


import argparse
import gzip
import os
import sys
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensembl.utils.archive import open_gz_file


def get_args():
    parser = argparse.ArgumentParser()
    # in files
    parser.add_argument(
        "fasta_dna",
        metavar="input.fna[.gz]",
        type=str,
        help="fasta (possibly gzipped) file with dna sequences to split",
    )
    # out
    parser.add_argument(
        "--out", metavar="chunks.fna", required=False, default=None, type=str, help="chunks output [STDOUT]"
    )
    parser.add_argument(
        "--individual_out_dir",
        required=False,
        default=None,
        type=str,
        help="output directory for writing files with individual chunks to (optional), "
        "if provided,`--out` value used as a filename prefix",
    )
    parser.add_argument(
        "--agp_out",
        metavar="chunks_contigs.agp",
        required=False,
        type=argparse.FileType("w", encoding="UTF-8"),
        help="AGP file with chunks to contigs mapping (optional)",
    )
    # meta_defaults
    parser.add_argument(
        "--chunk_size",
        metavar="100_000_000",
        required=False,
        type=int,
        default=100_000_000,
        help="maximum chunk size (should be greater then 50k)",
    )
    parser.add_argument(
        "--chunk_sfx",
        metavar="ens_chunk",
        required=False,
        type=str,
        default="ens_chunk",
        help="added to contig id before chunk number",
    )
    parser.add_argument(
        "--chunk_tolerance",
        metavar="ens_chunk",
        required=False,
        type=int,
        default=0,
        help="chunk size tolerance percentage (int). If the to-be-written chunk is longer "
        "than the defined chunk size (--chunk_size) by less than the specified "
        "tolerance percentage, it will not be split.",
    )
    parser.add_argument(
        "--n_seq",
        required=False,
        default=0,
        type=int,
        help="The fasta file will be split into chunks at positions of at least n N sequences. n number is specified here.",
    )
    parser.add_argument("--add_offset", required=False, action="store_true", help="0 based offset")

    args = parser.parse_args()
    if False and args.chunk_size < 50_000:
        parser.error(
            f"wrong '--chunk_size' value: '{args.chunk_size}'. should be greater then 50_000. exiting..."
        )
    if args.chunk_tolerance < 0:
        parser.error(
            f"wrong '--chunk_tolerance' value: '{args.chunk_tolerance}'. can't be less then 0. exiting..."
        )
    return args


def split_by_n(seq, split_re):
    """Split a string into chunks at the positions of N sequences.
    (return [ len(seq) ] if no split_re)"""
    seq_len = len(seq)
    if not split_re:
        return [seq_len]
    split_points = [m.end() for m in split_re.finditer(seq)]
    if not split_points or split_points[-1] != seq_len:
        split_points.append(seq_len)
    return split_points


def split_by_chunk_size(ends, chunk_size, tolerated_chunk_len=None):
    """Split array of end coordinates, to form chunks not longer then chunk_size (tolerated_chunk_len if speciefied)"""
    if tolerated_chunk_len is None or tolerated_chunk_len < chunk_size:
        tolerated_chunk_len = chunk_size
    result = []
    offset = 0
    for chunk_end in ends:
        chunk_len = chunk_end - offset
        if chunk_len > tolerated_chunk_len:
            # exclude starts, as they are 0 or pushed as previous chunk_ends
            result += list(range(offset, chunk_end, chunk_size))[1:]
        result.append(chunk_end)
        offset = chunk_end
    return result


## MAIN ##
def main():
    args = get_args()

    tolerated_chunk_len = args.chunk_size
    tolerated_chunk_len += args.chunk_size * args.chunk_tolerance // 100

    # make sure not used for args.n_seq <= 0
    n_split_regex = None
    if args.n_seq > 0:
        pattern = f"(N{{{args.n_seq},}})"
        n_split_regex = re.compile(pattern)

    # output_file to write sequences to
    file_prefix = ""
    if args.individual_out_dir:
        if not args.out:
            args.out = os.path.basename(args.fasta_dna)
        os.makedirs(args.individual_out_dir, exist_ok=True)
        file_prefix = os.path.join(args.individual_out_dir, args.out)
    else:
        out_file = open(args.out, "wt")

    # process input fasta
    fasta_file = args.fasta_dna
    with open_gz_file(fasta_file) as fasta:
        agp_lines = []
        print(
            f"spliting sequences from '{fasta_file}', chunk size {args.chunk_size:_}, splitting on {args.n_seq} Ns (0 -- disabled)",
            file=sys.stderr,
        )

        fasta_parser = SeqIO.parse(fasta, "fasta")
        for rec_count, rec in enumerate(fasta_parser, start=1):
            rec_name = str(rec.name)
            rec_len = len(rec)

            ends = split_by_n(str(rec.seq), n_split_regex)  # returns [ len(req.seq) ] if no regexp
            ends = split_by_chunk_size(ends, args.chunk_size, tolerated_chunk_len)

            offset = 0
            for chunk, chunk_end in enumerate(ends, start=1):
                chunk_name = f"{rec_name}_{args.chunk_sfx}_{chunk:03d}"
                chunk_file_name = f"{file_prefix}.{rec_count:03d}.{chunk:03d}.fa"
                if args.add_offset:
                    chunk_name += f"_off_{offset}"

                rec_from = offset + 1
                rec_to = chunk_end
                chunk_len = chunk_end - offset

                # form agp lines
                agp_line = f"{rec_name}\t{rec_from}\t{rec_to}\t{chunk}\tW\t{chunk_name}\t1\t{chunk_len}\t+"
                agp_lines.append(agp_line)

                # use agp lines as fasta description
                agp_line = agp_line.replace("\t", " ")
                print(f"dumping {chunk_name} AGP {agp_line}", file=sys.stderr)

                # get slice and put it out
                tmp_record = SeqRecord(
                    Seq(rec.seq[offset:chunk_end]), id=chunk_name, description=f"AGP {agp_line}", name=""
                )

                # if args.individual_files:
                #   ...
                if args.individual_out_dir:
                    with open(chunk_file_name, "wt") as out_file:
                        out_file.write(tmp_record.format("fasta"))
                else:
                    out_file.write(tmp_record.format("fasta"))

                del tmp_record
                #
                offset = chunk_end

        # dump AGP
        if args.agp_out:
            print("\n".join(agp_lines), file=args.agp_out)

        if not args.individual_out_dir:
            out_file.close()
    # main end


if __name__ == "__main__":
    main()
