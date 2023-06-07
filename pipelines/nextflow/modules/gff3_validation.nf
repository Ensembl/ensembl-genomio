// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

nextflow.enable.dsl=2

process GFF3_VALIDATION {

  //beforeScript 'module load libffi-3.3-gcc-9.3.0-cgokng6'
  tag "${gene_models}"
  label 'default'
  container "biocontainers/genometools:v1.5.10ds-3-deb_cv1"

  input:
    tuple val(gca), path (gene_models)

  output:
    tuple val(gca), path("*"), emit: gca

  script:
  """
  mv ${gene_models} gene_models.gff3.tmp 
  gt gff3 -tidy -sort -retainids -force -o ${gene_models} gene_models.gff3.tmp 
  gt gff3validator ${gene_models}
  """
}
