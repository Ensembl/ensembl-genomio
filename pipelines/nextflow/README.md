# Nextflow pipelines

## Installation
Please, refer to the "Installation" section of the [Nextflow pipelines document](/docs/nextflow.md#installation).

## Running
Please, refer to the "Running a pipeline" section of the [Nextflow pipelines document](/docs/nextflow.md#running_a_pipeline).
And to individual README(s) in '/pipelines/nextflow/workflows/' subfolders.

## Strange stuff with NextFlow pipelines
Please, refer to the "Strange things we met" section of the [Nextflow pipelines document](/docs/nextflow.md#Strange_things_we_met).

## Tree structure
This is an example of how the tree folder for Nextflow pipelines has been agreed to look like:
```
└── pipelines
    └── nextflow
        ├── conf
        │   └── base.config
        ├── modules
        │   └── fasta
        │   |   └── dump_fasta.nf
        |   └── utils.nf
        ├── subworkflows
        |    └── genome_validation
        │       ├── main.nf
        │       └── meta.yaml
        └── workflows
            ├── genome_loader
            │   ├── main.nf
            │   └── nextflow.config
            └── nextflow.config
```

## Pipelines test
To run the pipelines tests, run each test file with --git-aware:

```
pytest --git-aware pipelines/nextflow/tests/workflows/test_addition_prepare.yml
```

```
pytest --git-aware pipelines/nextflow/tests/workflows/test_genome_prepare.yml
```
