#!/usr/bin/env nextflow

process MANIFEST_STATS {
	tag "manifest_stats"                  
	label 'default'                

input:
	path manifest_dir
	val datasets

output:
	path manifest_dir

script:
 """
 manifest_stats --manifest_dir "$manifest_dir" --datasets_bin "$datasets"
 """
 }
