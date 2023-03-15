This is an example of how the tree folder for Nextflow pipelines has been agreed to look like:
```
└── pipelines
    └── nextflow
        ├── conf
        │   └── base.config
        ├── modules
        │   └── process_fasta.nf
        ├── subworkflows
        │   └── genome_validation.nf
        └── workflows
            ├── genome_loader
            │   ├── genome_loader.nf
            │   └── nextflow.config
            └── nexflow.config
```