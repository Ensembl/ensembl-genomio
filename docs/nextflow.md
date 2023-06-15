# Nextflow related documentation

## Installation

Define `NXF_HOME` env to use a nextlow home location instead of the default one (`$HOME/.nextflow`).
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
