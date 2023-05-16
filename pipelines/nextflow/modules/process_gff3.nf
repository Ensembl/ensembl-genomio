#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.merge_split_genes="False"

process PROCESS_GFF3 {
    tag "$gff3 - $task.attempt"
    label 'adaptive'

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

process gff3_validation {

  //beforeScript 'module load libffi-3.3-gcc-9.3.0-cgokng6'
  label 'default'
  container "biocontainers/genometools:v1.5.10ds-3-deb_cv1"

  input:
    path gene_models

  output:
    path "*"

  script:
  """
  mv ${gene_models} gene_models.gff3.tmp 
  gt gff3 -tidy -sort -retainids -force -o ${gene_models} gene_models.gff3.tmp 
  gt gff3validator ${gene_models}
  """
}


