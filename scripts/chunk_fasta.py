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
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_args():
  parser = argparse.ArgumentParser()
  # in files
  parser.add_argument("fasta_dna", metavar="input.fna[.gz]", type=str,
                      help="fasta (possibly gzipped) file with dna sequences to split")
  # out
  parser.add_argument("--out", metavar="chunks.fna", required = False,
                      type=argparse.FileType('w',  encoding='UTF-8'), default=sys.stdout,
                      help="chunks output [STDOUT]")
  parser.add_argument("--agp_out", metavar="chunks_contigs.agp", required = False,
                      type=argparse.FileType('w',  encoding='UTF-8'),
                      help="AGP file with chunks to contigs mapping (optional)")
  # meta_defaults
  parser.add_argument("--chunk_size", metavar="100_000_000", required = False,
                      type=int, default = 100_000_000, help="maximum chunk size (should be greater then 50k)")
  parser.add_argument("--chunk_sfx", metavar="ens_chunk", required = False,
                      type=str, default = "ens_chunk", help="added to contig id before chunk number")
  parser.add_argument("--chunk_tolerance", metavar="ens_chunk", required=False,
                      type=int, default=0, help="chunk size tolerance percentage. If the to-be-written chunk is longer "
                                                "than the defined chunk size (--chunk_size) by less than the specified "
                                                "tolerance percentage, it will not be split.")
  parser.add_argument("--n_seq", required = False, default = 0,
                      type=int, help="The fasta file will be split into chunks at positions of at least n N sequences. n number is specified here.")
  parser.add_argument("--add_offset", required = False, action="store_true",
                      help="0 based offset")

  args = parser.parse_args()
  if args.chunk_size < 50_000:
    parser.error(f"wrong '--chunk_size' value: '{args.chunk_size}'. should be greater then 50_000. exiting...")
  return args

def flatten(t):
    return [item for sublist in t for item in sublist]

def split_by_n(seq, n):
    """Split a string into chunks at the positions of N sequences"""
    pattern = f"(N{{{n},}})"
    regex = re.compile(pattern)
    split_points = [m.end() for m in regex.finditer(seq)]
    split_points.insert(0, 0)
    split_points.append(len(seq))
    return [(split_points[i], split_points[i+1]) for i in range(len(split_points)-1)]

def split_record_by_chunksize(starting_position, rec_len, chunk_size, tolerance):
    """Split a a range defined by a starting position and a range length (starting_position, rec_len)
    into chunks by a defined chunksize"""
    chunks = rec_len // chunk_size + (rec_len % chunk_size > 0)
    rec_from, rec_to, chunk_len = starting_position, starting_position+chunk_size, chunk_size
    chunks_list = []
    for chunk in range(1, chunks + 1):
        if (rec_to - starting_position) >= rec_len * (1 + tolerance/100):
            rec_to = starting_position + rec_len
        chunks_list.append((rec_from, rec_to))
        rec_from += chunk_len
        rec_to += chunk_len
    print(chunks_list)
    return chunks_list

## MAIN ##
def main():
  args = get_args()

  fasta_file = args.fasta_dna
  _open = fasta_file.endswith(".gz") and gzip.open or open
  with _open(fasta_file, 'rt') as fasta:
    agp_lines = []
    if not args.n_seq:
        print(f"spliting sequences from '{fasta_file}', chunk size {args.chunk_size:_}", file=sys.stderr)
    fasta_parser = SeqIO.parse(fasta, "fasta")
    for rec in fasta_parser:
      rec_len = len(rec)
      rec_name = str(rec.name)
      if args.n_seq != 0:
        n_chunks_list = split_by_n(str(rec.seq), int(args.n_seq))
        chunks_list = flatten([split_record_by_chunksize(x[0],(x[1]-x[0]),args.chunk_size, args.chunk_tolerance) for x in n_chunks_list])
        print(f"spliting {rec_name} ({rec_len:_} bp) into {len(chunks_list)} chunks on N sequences and specified chunk size.", file=sys.stderr)
      else:
        chunks_list = split_record_by_chunksize(0, rec_len, args.chunk_size, args.chunk_tolerance)
        print(f"spliting {rec_name} ({rec_len:_} bp) into {len(chunks_list)} chunks", file=sys.stderr)
      for chunk, chunk_set in enumerate(chunks_list):
        chunk+=1
        chunk_len = chunk_set[1] - chunk_set[0]
        rec_from = chunk_set[0] + 1
        rec_to = chunk_set[1]
        chunk_name = f"{rec_name}_{args.chunk_sfx}_{chunk:03d}"
        if (args.add_offset): chunk_name += f"_off_{rec_from - 1}"
        agp_line = f"{rec_name}\t{rec_from}\t{rec_to}\t{chunk}\tW\t{chunk_name}\t1\t{chunk_len}\t+"
        agp_lines.append(agp_line)
        agp_line = agp_line.replace("\t", " ")
        print(f"dumping {chunk_name} AGP {agp_line}", file = sys.stderr)
        tmp_record = SeqRecord(Seq(rec.seq[rec_from-1:rec_to]), id=chunk_name, description=f"AGP {agp_line}", name="")
        args.out.write(tmp_record.format("fasta"))
        del tmp_record
        #
    if args.agp_out:
      print("\n".join(agp_lines), file = args.agp_out)
  # main end


if __name__ == "__main__":
    main()
