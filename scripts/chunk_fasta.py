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
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
        "--out",
        metavar="chunks.fna",
        required=False,
        type=argparse.FileType("w", encoding="UTF-8"),
        default=sys.stdout,
        help="chunks output [STDOUT]",
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
    #
    args = parser.parse_args()
    if args.chunk_size < 50_000:
        parser.error(
            f"wrong '--chunk_size' value: '{args.chunk_size}'. should be greater then 50_000. exiting..."
        )
    return args


## MAIN ##
def main():
    args = get_args()

    fasta_file = args.fasta_dna
    _open = fasta_file.endswith(".gz") and gzip.open or open
    with _open(fasta_file, "rt") as fasta:
        agp_lines = []
        print(f"spliting sequences from '{fasta_file}', chunk size {args.chunk_size:_}", file=sys.stderr)
        fasta_parser = SeqIO.parse(fasta, "fasta")
        for rec in fasta_parser:
            rec_len = len(rec)
            rec_name = str(rec.name)
            chunks = rec_len // args.chunk_size + (rec_len % args.chunk_size > 0)
            print(f"spliting {rec_name} ({rec_len:_} bp) into {chunks} chunks", file=sys.stderr)
            #
            rec_from, rec_to, chunk_len = 1, args.chunk_size, args.chunk_size
            for chunk in range(1, chunks + 1):
                if rec_to > rec_len:
                    rec_to = rec_len
                    chunk_len = rec_len - rec_from + 1
                chunk_name = f"{rec_name}_{args.chunk_sfx}_{chunk:03d}"
                agp_line = f"{rec_name}\t{rec_from}\t{rec_to}\t{chunk}\tW\t{chunk_name}\t1\t{chunk_len}\t+"
                agp_lines.append(agp_line)
                #
                agp_line = agp_line.replace("\t", " ")
                print(f"dumping {chunk_name} AGP {agp_line}", file=sys.stderr)
                tmp_record = SeqRecord(
                    Seq(rec.seq[rec_from - 1 : rec_to]), id=chunk_name, description=f"AGP {agp_line}", name=""
                )
                args.out.write(tmp_record.format("fasta"))
                del tmp_record
                #
                rec_from += chunk_len
                rec_to += chunk_len
        if args.agp_out:
            print("\n".join(agp_lines), file=args.agp_out)
    # main end


if __name__ == "__main__":
    main()
