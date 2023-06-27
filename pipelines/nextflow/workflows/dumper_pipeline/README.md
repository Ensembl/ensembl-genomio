# A dumper pipeline
Creates seq\_regions.json and mysql dumps of the core database(s) from the given DBA host
with the given prefix (`--prefix`)
and/or with names matching regular expression (`--dbname_re`).
Output is created in the directory defined with `--output_dir` option (`./dumper_output` by default).

Sample output structure:
```
dumper_output/
├── coredb
│   └── EnsemblMetazoa
│       └── drosophila_melanogaster.sql.gz
└── metadata
    └── EnsemblMetazoa
        └── drosophila_melanogaster
            ├── manifest.json
            └── seq_region.json
```

N.B. `species.production_name`(s) are used to name output files.
Thus if you have several dbs using the same `species.production_name` output is incosistent.

Use `--help` to get more details.

See notes on [instalation](/docs/nextflow.md#Installation).

### Example invocation
```
CMD=<dba_alias>

mkdir -p data
pushd data
  data_dir=$(pwd)
  nextflow run \
    -w ${data_dir}/nextflow_work \
    ${ENSEMBL_ROOT_DIR}/ensembl-genomio/pipelines/nextflow/workflows/dumper_pipeline/main.nf \
    -profile lsf \
    $(${CMD} details script) \
    --dbname_re '^drosophila_melanogaster_\w+_57_.*$' \
    --output_dir ${data_dir}/dumper_output
popd
```

