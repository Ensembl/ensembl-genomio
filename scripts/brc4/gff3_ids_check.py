#!env python3

import re, sys, argparse

from BCBio import GFF

def check_gff_ids(gff_path):
    
    limits = {
            'gff_type' : ['gene'],
            }
    
    seqs_count = 0
    feats_count = 0
    
    prefixes = dict()
    for seq in GFF.parse(gff_path, limit_info=limits):
        seqs_count += 1
        for f in seq.features:
            feats_count += 1
            
            if f.type == 'gene':
                prefix = f.id.replace("gene-", "")
                prefix = prefix[0:4]
                
                if prefix in prefixes:
                    prefixes[prefix] += 1
                else:
                    prefixes[prefix] = 1
    
    print("GFF has %d seq_regions" % seqs_count)
    print("GFF has %d features" % feats_count)
    print("IDs prefixes:")
    for prefix in prefixes:
        print("%s : %d" % (prefix, prefixes[prefix]))

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Check stable ids in a gff3 file')
    
    parser.add_argument('--gff', type=str, required=True,
                help='GFF file to check')

    args = parser.parse_args()
    
    check_gff_ids(args.gff)

if __name__ == "__main__":
    main()
