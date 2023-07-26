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
import sys

from genmetaconf.metaconf import MetaConf
from genmetaconf.manifest import Manifest


def get_args():
    parser = argparse.ArgumentParser()
    # in files
    parser.add_argument(
        "--raw_meta_conf",
        metavar="meta/species",
        required=False,
        type=argparse.FileType("rt", encoding="UTF-8"),
        help="raw meta file with default settings",
    )
    parser.add_argument(
        "--gbff_file",
        metavar="species.gbff",
        required=False,
        type=str,
        help="gbff file to get species name, accession, taxid from",
    )
    parser.add_argument(
        "--fasta_dna",
        metavar="species_dna.fna",
        required=True,
        type=str,
        help="fasta file with dna sequences",
    )
    parser.add_argument(
        "--fasta_pep",
        metavar="species_pep.fa",
        required=False,
        type=str,
        help="fasta file with protein sequences",
    )
    parser.add_argument(
        "--gff_file", metavar="species.gff3", required=False, type=str, help="gff file with  gene models"
    )
    parser.add_argument(
        "--fann_file",
        metavar="functional_annotation.json",
        required=False,
        type=str,
        help="json file with functional annotation",
    )
    parser.add_argument(
        "--seq_region_raw",
        metavar="seq_region_raw.json",
        required=False,
        type=str,
        help="seq_region raw json file to patch",
    )
    parser.add_argument(
        "--seq_region_genbank",
        metavar="seq_region_genbank.json",
        required=False,
        type=str,
        help="ad-hoc seq_region raw json file based on the genbank gff3",
    )
    parser.add_argument(
        "--seq_region_syns",
        metavar="seq_region_syns.tsv",
        required=False,
        type=str,
        help="ad-hoc seq_region synomyms [source, optionally] tab-separated file (known_name \\t syn [ \\t source ])",
    )
    parser.add_argument(
        "--asm_rep_file",
        metavar="species_assembly_report.txt",
        required=False,
        type=str,
        help="GenBank assembly report to get seq_region syns from",
    )
    # out
    parser.add_argument(
        "--meta_out",
        metavar="data/metadata/species",
        required=False,
        type=argparse.FileType("w", encoding="UTF-8"),
        default=sys.stdout,
        help="stats output [STDOUT]",
    )
    parser.add_argument(
        "--data_out_dir", metavar="data/metadata", required=True, type=str, help="dir to store files into"
    )
    parser.add_argument(
        "--genome_conf", metavar="genome.json", required=True, type=str, help="genome json file output"
    )
    parser.add_argument(
        "--seq_region_conf",
        metavar="seq_region.json",
        required=False,
        type=str,
        help="seq_region json file output",
    )
    parser.add_argument(
        "--manifest_out", metavar="manifest.json", required=True, type=str, help="manifest file output"
    )
    # meta_defaults
    parser.add_argument(
        "--assembly_version",
        metavar="1",
        required=False,
        type=int,
        default=1,
        help="assembly.version default",
    )
    parser.add_argument(
        "--species_division",
        metavar="EnsemblMetazoa",
        required=False,
        type=str,
        default="EnsemblMetazoa",
        help="species.division default",
    )
    parser.add_argument(
        "--genebuild_method",
        metavar="import",
        required=False,
        type=str,
        default="import",
        help="genebuild.method default",
    )
    parser.add_argument(
        "--genebuild_method_display",
        metavar="Import",
        required=False,
        type=str,
        default="Import",
        help="genebuild.method_display default",
    )
    parser.add_argument(
        "--genebuild_level",
        metavar="toplevel",
        required=False,
        type=str,
        default="toplevel",
        help="genebuild.level default",
    )
    parser.add_argument(
        "--syns_src",
        metavar="GenBank",
        required=False,
        type=str,
        default="GenBank",
        help="seq region syns source default",
    )
    parser.add_argument(
        "--default_genetic_code",
        metavar="1",
        required=False,
        type=int,
        default=1,
        help="default genetic code variant",
    )
    parser.add_argument(
        "--default_circular",
        action="store_true",
        required=False,
        default=False,
        help="assume contigs are circular by default",
    )
    parser.add_argument(
        "--generate_species_aliases",
        action="store_true",
        required=False,
        default=False,
        help="add generated species aliases",
    )
    #
    args = parser.parse_args()
    return args


## MAIN ##
def main():
    args = get_args()

    meta = MetaConf(args.raw_meta_conf, args.generate_species_aliases)
    meta.merge_from_asm_rep(args.asm_rep_file)
    meta.merge_from_gbff(args.gbff_file)
    meta.update_derived_data(
        {
            "assembly.version": args.assembly_version,
            "genebuild.method": args.genebuild_method,
            "genebuild.method_display": args.genebuild_method_display,
            "genebuild.level": args.genebuild_level,
            "species.division": args.species_division,
        },
        update_annotation_related=(args.gff_file is not None),
    )
    meta.dump(args.meta_out)

    meta.dump_genome_conf(args.genome_conf)
    meta.dump_seq_region_conf(
        args.seq_region_conf,
        fasta_file=args.fasta_dna,
        asm_rep_file=args.asm_rep_file,
        seq_region_raw=args.seq_region_raw,
        seq_region_genbank=args.seq_region_genbank,
        seq_region_syns=args.seq_region_syns,
        syns_src=args.syns_src,
        default_genetic_code=args.default_genetic_code,
        default_circular=args.default_circular,
    )
    # gen manifest
    manifest = Manifest(
        {
            "fasta_dna": args.fasta_dna,
            "fasta_pep": args.fasta_pep,
            "genome": args.genome_conf,
            "seq_region": args.seq_region_conf,
            "functional_annotation": args.fann_file,
            "gff3": args.gff_file,
        }
    )
    manifest.dump(args.manifest_out, outdir=args.data_out_dir)


if __name__ == "__main__":
    main()
