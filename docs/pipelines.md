# Ensembl Genomio Pipelines:

## Genomio prepare pipeline
_Module [Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_prepare_conf]_

**Genome prepare pipeline for BRC/Metazoa**

#### Description
Retrieve data for a genome from INSDC and prepare the following files in a separate folder
for each genome:

- FASTA for DNA sequences
- FASTA for protein sequences
- GFF gene models
- JSON functional annotation
- JSON seq_region
- JSON genome
- JSON manifest

The JSON files follow the schemas defined in the `src/python/ensembl/io/genomio/data/schemas` folder.

These files can then be fed to the Genome loader pipeline.

### How to run

```
init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_prepare_conf \
    --host $HOST --port $PORT --user $USER --pass $PASS \
    --hive_force_init 1 \
    --pipeline_dir temp/prepare \
    --data_dir $INPUT \
    --output_dir $OUTPUT \
    ${OTHER_OPTIONS}
```

### Parameters

| option | default value |  meaning |
| - | - | - |
| `--pipeline_name` | brc4_genome_prepare | name of the hive pipeline
| `--pipeline_dir` | | temp directory for this pipeline run
| `--data_dir` | | directory with json files for each genome to prepare, following the format set by `src/python/ensembl/io/genomio/data/schemas/genome.json`
| `--output_dir` | | directory where the prepared files are to be stored
| `--merge_split_genes` | 0 | Sometimes the gene features are split in a gff file. Ensembl expects genes to be contiguous, so this option merge the parts into 1.
| `--exclude_seq_regions` |  | Do not include those seq_regions (apply to all genomes, this should be seldom used)
| `--validate_gene_id` | 0 | Enforce a strong gene ID pattern (replace by GeneID if available)
| `--ensembl_mode` |  0 | By default, set additional metadata for BRC genomes. With this parameter, use vanilla Ensembl metadata.
