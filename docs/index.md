# Ensembl GenomIO

A repository dedicated to pipelines used to turn basic genomic data into formatted 
Ensembl core databases. Also allow users to dump core databases into various formats.

File formats handled : FastA, GFF3, JSON (*following BRC4 specifications*).

## Contents

Check out [installation](install.md) section for further information on how 
to install the project.

1. [Install](install.md)
2. [Usage](usage.md)
3. [Code of Conduct](code_of_conduct.md)
4. [Code reference](reference/ensembl/io/genomio)

## Ehive pipelines

Check out the [usage](usage.md) section for further information of requirements to
run ensembl-genomio pipelines.

1. __Genome loader__: Creates an Ensembl core database from a set of flat files.
2. __Genome dumper__: Dumps flat files from an Ensembl core database.

## Nextflow pipelines

1. __Additional seq prepare__: BRC/Ensembl metazoa pipeline. Preparation of genome data loading files for new sequence(s) to existing species databases.  
2. __Genome Prepare__: BRC/Ensembl metazoa pipeline. Retrieve data for genome(s), obtained from INSDC and RefSeq, validate and prepare GFF3, FASTA, JSON files for each genome accession.

## License

Software as part of [Ensembl GenomIO](https://github.com/Ensembl/ensembl-genomio) is distributed under the [Apache-2.0 License](https://www.apache.org/licenses/LICENSE-2.0.txt).
