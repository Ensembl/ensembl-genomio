#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process extract_from_gb {
    tag "$gb_file"
    label 'default'

    input:
    val prefix
    val production_name
    path gb_file 

    output:
    path "*.gff", emit: gene_gff
    path "genome.json", emit: genome
    path "seq_region.json", emit: seq_regions
    path "dna.fasta", emit: dna_fasta
    path "pep.fasta", emit: pep_fasta   //is this optional?

    script:
    """
    extract_from_gb --prefix ${prefix} --prod_name ${production_name} --gb_file ${gb_file}
    """
}

