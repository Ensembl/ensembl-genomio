# Nextflow related documentation

## Installation
If you don't have an installed environment or you don't have nextflow itself, here's one of the ways to install it.

Define [`NXF_HOME` env variable](https://www.nextflow.io/docs/latest/config.html#environment-variables) to use a nextlow home location instead of the default one (`$HOME/.nextflow`).
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

Configure the environment you're using if you haven't done so yet.
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
Or use `nexflow -e.NXF_WORK=...` approach.
Ideally, should be overridable by the `-work-dir` (`-w`) option of `nextflow run`

## Running a pipeline

