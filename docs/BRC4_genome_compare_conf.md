# **Compare genomes pipeline**
*MODULE*: *[Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_compare_conf](https://github.com/Ensembl/ensembl-genomio/blob/main/src/perl/Bio/EnsEMBL/Pipeline/PipeConfig/BRC4_genome_compare_conf.pm)*

## **Overview**
-----
This pipeline is used for a sequence-level comparison of an assembly with INSDC and provides a detailed report on the discrepancies. The following steps are performed:

  1. Download the files for the corresponding assembly from INSDC
  2. Retrieve metadata seq.json and fasta files from the database
  3. Compare the fasta files
       - compare the sequence ids
       - compare the sequence 
       - MD5 check of the files
       - identify organellar sequences in both assemblies
  4. Report the results

## **Prerequisites**
A [registry file](https://www.ensembl.org/info/docs/api/registry.html) to connect to the database.

## **Example of how to create the pipeline**
```
init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_compare_conf \
  --host $HOST --port $PORT --user $USER --pass $PASS \
  --registry $REGISTRY \
  --hive_force_init 1 \
  --output_dir $OUTPUT_DIR \
  --tmp_dir temp/compare \
  --species acanthamoeba_astronyxis_gca000826245
```

### **PARAMETERS**
| Options | Type | Default value | Mandatory | Description |
| - | - |  - |  - | - |
| `--registry` | file |  | yes | service that connects to the database |
| `--pipeline_name` | str | brc4_genome_compare |   optional| name of the hive pipeline |
| `--hive_force_init` | int |  | yes | drop and create the hive pipeline from scratch | 
| `--output_dir`      | dir |   ./output                     | optional| directory to store the result | 
| `--tmp_dir`         | dir |   ./tmp                     | optional| temp directory for downloaded files |
| `--species`         | str |                        | yes| species (one or multiple) to process (production name) |
| `--run_all`         | int |     0                   | yes| process all the species in the registry | 
| `--email`           | str |  $USER.ebi.ac.uk                     | optional| a summary is emailed when the pipeline is complete | 

Note:
Either use ```--species``` to run one or multiple species separately or ```--run_all 1``` for all the species in the database.
Currently this pipeline is only used to compare with Genbank assembly. 
## **RESULT**
---------------------------------------------------
Generates 3 files:
  - report.log: A tab-delimited file containing a summary of the compared sequences between the INSDC assembly/assemblies and the database(s)
    - The report contains 13 columns: 
      - species: name of the species 
      - accession: GCA accession
      - seq_count_1: total number of sequences in INSDC
      - seq_count_2: total number of sequences in the database
      - num_diff_seq: the total number of sequences that differ between INSDC and the database
      - common: the total number of common sequences between INSDC and the database
      - only1: count of sequences found only in INSDC
      - only2: count of sequences found only in the database
      - max_only1: a total of the sequence length in only1
      - max_only2: a total of the sequence length in only2
      - other_locations: the total count of organellar genomes 
      - summary (mismatch or identical)
      - organellar_summary
  
  - species_fasta_dna.map: A [JSON schema](https://github.com/Ensembl/ensembl-genomio/blob/main/src/python/ensembl/io/genomio/data/schemas/seq_region.json) file containing metadata of the common sequences
  - species_fasta_dna.log: Detailed report of mismatched sequences
