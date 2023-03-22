#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process download_genbank {
    label "Sequence_genbank_file"
    tag "${accession}"
    label 'default'

    input:
        val accession

    output:
       path "*.gb"

    script:
    """
    download_genbank --accession ${accession}
    """
}

process extract_from_gb {
    tag "$gb_file"
    label 'default'

    input:
    val prefix
    val prod_name
    path gb_file 

    output:
    path "*.gff", emit: gene_gff
    path "genome.json", emit: genome
    path "seq_region.json", emit: seq_regions
    path "dna.fasta", emit: dna_fasta
    path "pep.fasta", emit: pep_fasta   //is this optional?

    script:
    """
    extract_from_gb --prefix ${prefix} --prod_name ${prod_name} --gb_file ${gb_file}
    """
}

