# Nextflow related documentation

## Installation
If you do not have an installed environment or you don't have nextflow itself, here is one of the ways to install it.

Define [`NXF_HOME` env variable](https://www.nextflow.io/docs/latest/config.html#environment-variables) to use a nextflow home location instead of the default one (`$HOME/.nextflow`).
Everything else is unchanged from the default Nextflow installation instructions on [https://www.nextflow.io/index.html#GetStarted](https://www.nextflow.io/index.html#GetStarted).

```
# add NXF_HOME env
export NXF_HOME=$(pwd)/dot.nextflow # or whatever

# get nextflow and install almost like here: https://www.nextflow.io/index.html#GetStarted
wget -O - https://get.nextflow.io  > nextflow.install.bash

# review and run
cat nextflow.install.bash | bash -i 2>&1 | tee nextflow.install.log

# run test, see https://www.nextflow.io/index.html#GetStarted
./nextflow run hello
```

Configure the environment you are using if you have not done so yet.
Don't forget to add `NXF_HOME`, patch `PATH` and export them.
```
# fix env variables, i.e.:
export NXF_HOME=$(pwd)/dot.nextflow
export PATH=$(pwd):$PATH
```

If you wish, you can set `NXF_WORK` env to be used by `nextflow`.
```
export NXF_WORK=...
```
Or use `nextflow -e.NXF_WORK=...` approach.
Ideally, should be overridable by the `-work-dir` (`-w`) option of `nextflow run`

## Running a pipeline
Once you have production (and nextflow) env ready, you can run pipelines.
I.e.
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

Try to invoke pipelines with `--help` option to get insight on how to run them.

## Strange things we met

### Channel is not forked, only one operation on stream is allowed
#### Symptoms:
When running a stage or a subworkflow on a channel with a single element we expect stream to be forked, allowing us to seed several task at a time.
```
// create that channel with a single element
//   calls read_json(...) in turn, see below
dbs = from_read_json(...)

DUMP_SQL(..., dbs, ...)
DUMP_METADATA(..., dbs, ...)
```

Instead pipeline dies with
```
Caused by: Cannot load from object array because "this.keys" is null
```
and when printing this object (`dbs` in this case, with `println "db: ${db}"`), we see it dict surrounded by the curly brackets like this
```
{..., "db_name":"some_db_name", ...}
```
instead of this (with square brackets)
```
[..., "db_name":"some_db_name", ...]
```

#### Reason / solution
In our case we used the `read_json` function similar to this one:
```
def read_json(json_path) {
    slurp = new JsonSlurper()
    json_file = file(json_path)
    text = json_file.text
    return slurp.parseText(text) // <-- problem here
}
```
that returned some kind of a lazy evaluator/iterator/whatever(not sure).

Replacing `return slurp.parseText(text)` with
```
    not_a_lazy_val = slurp.parseText(text)
    return not_a_lazy_val
```
did help.
