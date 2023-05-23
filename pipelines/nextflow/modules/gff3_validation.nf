#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GFF3_VALIDATION {

  //beforeScript 'module load libffi-3.3-gcc-9.3.0-cgokng6'
  tag "${gene_models}"
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
