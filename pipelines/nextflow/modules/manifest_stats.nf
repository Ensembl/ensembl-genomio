#!/usr/bin/env nextflow

process MANIFEST_STATS {
	tag "manifest_stats"                  
	label 'default'                
	publishDir "$manifest_dir", mode: 'copy' 

input:
	path manifest_dir
	val accession  
	val datasets

output:
	path "$manifest_dir/stats.txt"

script:
 """
 manifest_stats --manifest_dir $manifest_dir --accession $accession --datasets_bin $datasets 
 """
 }
 
