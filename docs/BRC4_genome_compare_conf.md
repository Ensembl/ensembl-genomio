# **Compare genomes pipeline**
*MODULE*: *[Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_compare_conf](https://github.com/Ensembl/ensembl-genomio/blob/main/lib/perl/Bio/EnsEMBL/Pipeline/PipeConfig/BRC4_genome_compare_conf.pm)*

## **Overview**
-----
This pipeline is used for a sequence-level comparison of an assembly with INSDC and provides a detailed report on the discrepencies. The following steps are performed:

  1. Download the files for the corresponding assembly from INSDC
  2. Retreive metadata seq.json and fasta files from the database
  3. Compare the fasta files
       - compare the sequence ids
       - compare the sequence 
       - MD5 check of the files
       - identify organellar sequences in both assemblies
  4. Report the results

## **How to run**:
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
| Options | type | default value | mandatory |  description |
| - | - |  - |  - | - |
| `--registry` | dir |  | yes | service that connects to the database |
| `--pipeline_name` | str | brc4_genome_compare |   optional| name of the hive pipeline |
| `--hive_force_init` | int |  | yes | drop and run the hive pipeline from scratch | 
| `--output_dir`      | dir |   ./output                     | optional| directory to store the result | 
| `--tmp_dir`         | dir |   ./tmp                     | optional| temp directory for dowloaded files |
| `--species`         | str |                        | yes| species to process (production name) |
| `--run_all`         | int |     0                   | yes| process all the species in the registry | 
| `--email`           | str |  $USER.ebi.ac.uk                     | optional| a summary is emailed when the pipeline is complete | 

Note:
Either use --species  to run one or multiple species separately or --run_all 1 for all the species in the database.
Currently this pipeline is only used to compare with Genbank assembly. 
## **RESULT**
---------------------------------------------------
Generates 3 files:
  - report.log: A tab-delimited file containing a summary of the compared sequences
  - species_fasta_dna.map: A JSON format file containing metadata of the common sequences
  - species_fasta_dna.log: Detailed report of mismatched sequences










