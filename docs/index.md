# [Ensembl-genomIO (v0.1)](https://github.com/Ensembl/ensembl-genomio)

*Ensembl-genomIO Base Library Documentation*

A repository dedicated to pipelines used to turn basic genomic data into formatted 
Ensembl core databases. Also allow users to dump core databases into various formats.

File formats handled : FastA, GFF3, JSON (*following BRC4 specifications*).

Contents
--------
Check out [installation](install.md) section for further information on how 
to install the project.

1. [Usage](usage.md)
2. [Install](install.md)
3. [License](license.md)

Ehive pipelines
-------------------------------------------
Check out the [usage](usage.md) section for further information of requirements to
run ensembl-genomio pipelines.

1. __Genome loader__: Creates an Ensembl core database from a set of flat files.
2. __Genome dumper__: Dumps flat files from an Ensembl core database.

Nextflow pipelines
-------------------------------------------
1. __Additional seq prepare__: BRC/Ensembl metazoa pipeline. Preparation of genome data loading files for new sequence(s) to existing species databases.  
2. __Genome Prepare__: BRC/Ensembl metazoa pipeline. Retrieve data for genome(s), obtained from INSDC and RefSeq, validate and prepare GFF3, FASTA, JSON files for each genome accession.


## Project layout
	src/ensembl/
	├── brc4
	│   └── runnable
	│       ├── compare_fasta.py
	│       ├── compare_report.py
	│       ├── core_server.py
	│       ├── download_genbank.py
	│       ├── dump_stable_ids.py
	│       ├── extract_from_gb.py
	│       ├── fill_metadata.py
	│       ├── gff3_specifier.py
	│       ├── integrity.py
	│       ├── json_schema_factory.py
	│       ├── load_sequence_data.py
	│       ├── manifest.py
	│       ├── manifest_stats.py
	│       ├── prepare_genome.py
	│       ├── read_json.py
	│       ├── say_accession.py
	│       └── seqregion_parser.py
	└── io
	    └── genomio
	        ├── assembly
	        │   └── download.py
	        ├── database
	        │   └── factory.py
	        ├── events
	        │   ├── dump.py
	        │   ├── format.py
	        │   └── load.py
	        ├── fasta
	        │   └── process.py
	        ├── genbank
	        │   ├── download.py
	        │   └── extract_data.py
	        ├── genome_metadata
	        │   ├── dump.py
	        │   ├── extend.py
	        │   └── prepare.py
	        ├── genome_stats
	        │   ├── compare.py
	        │   └── dump.py
	        ├── gff3
	        │   ├── extract_annotation.py
	        │   └── process.py
	        ├── manifest
	        │   ├── check_integrity.py
	        │   ├── compute_stats.py
	        │   └── generate.py
	        ├── schemas
	        │   └── json
	        │       ├── factory.py
	        │       └── validate.py
	        ├── seq_region
	        │   ├── dump.py
	        │   └── prepare.py
	        └── utils
	            ├── archive_utils.py
	            └── json_utils.py
