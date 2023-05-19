#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.merge_split_genes="False"

process PROCESS_GFF3 {
    tag "$gff3 - $task.attempt"
    label 'variable_2_8_32'

    input:
    path gff3
    path genome
    val gca

    output:
    path "*.gff3", emit: gene_models
    path "*.json", emit: functional_annotation
    val gca, emit: gca

    script:
    """
    process_gff3 --in_gff_path ${gff3} --genome_data ${genome}
    """
}
