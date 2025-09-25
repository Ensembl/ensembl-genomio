# Ensembl GenomIO

[![License](https://img.shields.io/badge/license-Apache_2.0-blue.svg)](https://github.com/Ensembl/ensembl-genomio/blob/main/LICENSE)
[![Coverage](https://plantazoa.gitdocs.ebi.ac.uk/ensembl-genomio/coverage-badge.svg)](https://plantazoa.gitdocs.ebi.ac.uk/ensembl-genomio/)
[![CI](https://img.shields.io/github/checks-status/Ensembl/ensembl-genomio/main?label=CI)](https://gitlab.ebi.ac.uk/plantazoa/ensembl-genomio/-/pipelines)
[![Release](https://img.shields.io/pypi/v/ensembl-genomio)](https://pypi.org/project/ensembl-genomio)

Pipelines to turn basic genomic data into Ensembl cores and back.

This is a multilanguage (Perl, Python) repo providing eHive pipelines and various scripts (see below) to prepare genomic data and load it as [Ensembl core database](http://www.ensembl.org/info/docs/api/core/index.html) or to dump such core databases as file bundles.

Bundles themselves consist of genomic data in various formats (e.g. fasta, gff3, json) and should follow the corresponding [specification](https://github.com/Ensembl/ensembl-genomio/blob/main/docs/BRC4_genome_loader.md#input-data).


## Installation and configuration

This repository is publicly available in [PyPI](https://pypi.org), so it can be easily installed with your preferred Python package manager, e.g.:

```bash
pip install ensembl-genomio
```

### Prerequisites

Pipelines are intended to be run inside the Ensembl production environment. Please, make sure you have all the proper credential, keys, etc. set up.

### Get repo and install

Clone:
```
git clone git@github.com:Ensembl/ensembl-genomio.git
```

Install the python part (of the pipelines) and test it:
```
pip install ./ensembl-genomio
# And test it has been installed correctly
python -c 'import ensembl.io.genomio'
```

Update your perl envs (if you need to)
```
export PERL5LIB=$(pwd)/ensembl-genomio/src/perl:$PERL5LIB
export PATH=$(pwd)/ensembl-genomio/scripts:$PATH
```

### Optional installation

If you need to install "editable" Python package use '-e' option
```
pip install -e ./ensembl-genomio
```

To install additional dependencies (e.g. `[docs]` or `[cicd]`) provide `[<tag>]` string, e.g.:
```
pip install -e ./ensembl-genomio[cicd]
```

For the list of tags see `[project.optional-dependencies]` in [pyproject.toml](https://github.com/Ensembl/ensembl-genomio/blob/main/pyproject.toml).

### Additional steps to use automated generation of the documentation

- Install python part with the `[docs]` tag
- Change into repo dir
- Run `mkdocs build` command

```
git clone git@github.com:Ensembl/ensembl-genomio.git
cd ./ensembl-genomio
pip install -e .[docs]
mkdocs build
```

###  Nextflow installation

Please, refer to the "Installation" section of the [Nextflow pipelines document](https://github.com/Ensembl/ensembl-genomio/blob/main/docs/nextflow.md#installation).

## Pipelines

### Initialising and running eHive-based pipelines

Pipelines are derived from [`Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf`](https://github.com/Ensembl/ensembl-hive/blob/version/2.6/modules/Bio/EnsEMBL/Hive/PipeConfig/HiveGeneric_conf.pm),
or from [`Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf`](https://github.com/Ensembl/ensembl-hive/blob/version/2.6/modules/Bio/EnsEMBL/Hive/PipeConfig/EnsemblGeneric_conf.pm),
of from [`Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf`](https://github.com/Ensembl/ensembl-production-imported/blob/main/src/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/EGGeneric_conf.pm) (see [documentation](https://github.com/Ensembl/ensembl-production-imported/blob/main/docs/EGGeneric.md)).

And the same perl class prefix used for every pipeline:
  `Bio::EnsEMBL::EGPipeline::PipeConfig::` .

N.B. Don't forget to specify `-reg_file` option for the `beekeeper.pl -url $url -reg_file $REG_FILE -loop` command.

```
init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_loader_conf
    $($CMD details script) \
    -hive_force_init 1\
    -queue_name $SPECIFIC_QUEUE_NAME \
    -registry $REG_FILE \
    -pipeline_tag "_${PIPELINE_RUN_TAG}" \
    -ensembl_root_dir ${ENSEMBL_ROOT_DIR} \
    -dbsrv_url $($CMD details url) \
    -proddb_url "$($PROD_SERVER details url)""$PROD_DBNAME" \
    -taxonomy_url "$($PROD_SERVER details url)""$TAXONOMY_DBNAME" \
    -release ${RELEASE_VERSION} \
    -data_dir ${INPUT_DIR}/manifests_dir/ \
    -pipeline_dir $OUT_DIR/loader_run \
    ${OTHER_OPTIONS} \
    2> $OUT_DIR/init.stderr \
    1> $OUT_DIR/init.stdout

SYNC_CMD=$(cat $OUT_DIR/init.stdout | grep -- -sync'$' | perl -pe 's/^\s*//; s/"//g')
# should get something like
#   beekeeper.pl -url $url -sync

LOOP_CMD=$(cat $OUT_DIR/init.stdout | grep -- -loop | perl -pe 's/^\s*//; s/\s*#.*$//; s/"//g')
# should get something like
#   beekeeper.pl -url $url -reg_file $REG_FILE -loop

$SYNC_CMD 2> $OUT_DIR/sync.stderr 1> $OUT_DIR/sync.stdout
$LOOP_CMD 2> $OUT_DIR/loop.stderr 1> $OUT_DIR/loop.stdout
```

### List of the pipelines

| Pipeline name | Description | Document | Comment | Module |
| - | - | - | - | - |
| BRC4_genome_loader | creates an [Ensembl core database](http://www.ensembl.org/info/docs/api/core/index.html) from a set of flat files or adds ad-hoc (i.e. organellas) sequences to the existing core  | [BRC4_genome_loader](https://github.com/Ensembl/ensembl-genomio/blob/main/docs/BRC4_genome_loader.md) | | [Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_loader_conf](https://github.com/Ensembl/ensembl-genomio/blob/main/src/perl/Bio/EnsEMBL/Pipeline/PipeConfig/BRC4_genome_loader_conf.pm)
| BRC4_genome_dumper | | | | |
| | | | | |
| BRC4_genome_prepare | | | | |
| BRC4_addition_prepare | | | | |
| BRC4_genome_compare | | | | |
| | | | | |
| LoadGFF3 | | | | |
| LoadGFF3Batch | | | | |


### Scripts

* [trf_split_run.bash](https://github.com/Ensembl/ensembl-genomio/blob/main/scripts/trf_split_run.bash) -- a trf wrapper with chunking support to be used with [ensembl-production-imported DNAFeatures pipeline](https://github.com/Ensembl/ensembl-production-imported/tree/main/src/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/DNAFeatures_conf.pm) (see [docs](https://github.com/Ensembl/ensembl-genomio/blob/main/docs/trf_split_run.md))

## CI/CD bits

As for now some [Gitlab CI](https://docs.gitlab.com/ee/ci/) pipelines introduced to keep things in shape. Though, this bit is in constant development. Some documentation can be found in [docs for GitLab CI/CD](https://github.com/Ensembl/ensembl-genomio/blob/main/docs/cicd_gitlab.md)

## Various docs

See [docs](https://github.com/Ensembl/ensembl-genomio/blob/main/docs)

## Unit testing

The Python part of the codebase has now unit tests available to test each module. Make sure you have installed this repository's `[cicd]` dependencies (via `pip install ensembl-genomio[cicd]`) before continuing.

Running all the tests in one go is as easy as running `pytest` **from the root of the repository**. If you also want to measure, collect and report the code coverage, you can do:
```bash
coverage run -m pytest
coverage report
```

You can also run specific tests by supplying the path to the specific test file/subfolder, e.g.:
```bash
pytest lib/python/tests/test_schema.py
```

## Acknowledgements

Some of this code and documentation is inherited from the [EnsemblGenomes](https://github.com/EnsemblGenomes) and other [Ensembl](https://github.com/Ensembl) projects. We appreciate the effort and time spent by developers of the [EnsemblGenomes](https://github.com/EnsemblGenomes) and [Ensembl](https://github.com/Ensembl) projects. 

Thank you!
