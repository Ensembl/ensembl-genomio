# How to run those pipelines

## Genome dumper
This pipeline dumps files from Ensembl core databases following the specifications for file transfers between EBI and EuPathDB for BRC4.

### Init command
init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_dumper_conf \
  --schema_dir .../ensembl-genomio/schemas \
  $($SERVER details hive) \
  --release 99 \
  --output_dir $OUTPUT \
  --tmp_dir $TMP \
  --registry $REGISTRY \
  --run_all 1

### Input
The pipeline requires access to Ensembl core databases listed in the registry, and a list of species to run (or run_all).

### Output
The output dir will contain one directory per species (named with their production_name), containing all files (fasta, gff, json) and a manifest for those files.


## Genome loader
TODO
