# Nextflow pipelines

## Installation
Please, refer to the "Installation" section of the [Nextflow pipelines document](/docs/nextflow.md#installation).

## Running
Please, refer to the "Running a pipeline" section of the [Nextflow pipelines document](/docs/nextflow.md#running_a_pipeline).
And to individual README(s) in '/pipelines/nextflow/workflows/' subfolders.

### Trace and log reports
Add `-with-trace` to your `nextflow run ...` command to activate the generation of a trace file with useful information such as resource usage of each job and their hash directory. More information can be found in the [Nextflow documentation](https://www.nextflow.io/docs/latest/tracing.html).

Most of our jobs in our pipelines have the `--verbose` flag set as default, providing a more detailed log. You can generate a single log report with the information from every job by running:
```bash
nextflow log <run name> -t log_template.md > log_report.md
```

However, sometimes this report may be too large to be useful, so the previous command includes a `-filter` option, e.g. `-filter 'name =~ /DOWNLOAD_ASM_DATA.*/'` will provide a log report only for those jobs whose name starts with `DOWNLOAD_ASM_DATA`.

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
