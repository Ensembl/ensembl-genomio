# ensembl-genomio
Pipelines to turn basic genomic data into Ensembl cores and back

This is a mulitlanguage (Perl, Python) repo providing eHive pipelines
and various scripts (see below) to prepare genomic data and load it as [Ensembl core database](http://www.ensembl.org/info/docs/api/core/index.html) or to dump such core databases as file bundles.

The files composing data bundle are of various formats (fasta, gff3, json) and should follow the corresponding [specification](docs/BRC4_genome_loader.md#input).


## Installation and configuration
This repo 

### Prerequisites
Pipelines are intended to be run inside the Ensembl production environment.
Please, make sure you have all the proper credential, keys, etc. set up.

### Get repo and install

Clone:
```
git clone git@github.com:Ensembl/ensembl-genomio.git 
```

Install the python part (of the pipelines) and test it:
```
pip install ./ensembl-genomio

# test
python -c 'import ensembl.brc4.runnable.read_json'
```

Update your perl envs (if you need to)
```
export PERL5LIB=$(pwd)/ensembl-genomio/lib/perl:$PERL5LIB
export PATH=$(pwd)/ensembl-genomio/scripts:$PATH
```

### Optional installation
If you need to install "editable" python package use '-e' option
```
pip install -e ./ensembl-genomio
```

To install additional dependencies (e.g. `[doc]` or `[dev]`) provide `[<tag>]` string. I.e.
```
pip install -e ./ensembl-genomio[dev]
```

For the list of tags see `[project.optional-dependencies]` in [pyproject.toml](./pyproject.toml). 


### Additional steps to use automated genertaion of the documentation (part of it)
Install python part with the `[doc]` or `[dev]` tag.
Change into repo dir
Run doc build script.

```
git clone git@github.com:Ensembl/ensembl-genomio.git 
pip install -e ./ensembl-genomio[doc]

cd ./ensembl-genomio

# build docs
./scripts/setup/docs/build_sphinx_docs.sh
```

## Pipelines

### Initialising and running eHive-based pipelines

Pipelines are derived from [`Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf`](https://github.com/Ensembl/ensembl-hive/blob/version/2.6/modules/Bio/EnsEMBL/Hive/PipeConfig/HiveGeneric_conf.pm),
or from [`'Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf`](https://github.com/Ensembl/ensembl-hive/blob/version/2.6/modules/Bio/EnsEMBL/Hive/PipeConfig/EnsemblGeneric_conf.pm),
of from [`Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf`](https://github.com/Ensembl/ensembl-production-imported/blob/main/lib/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/EGGeneric_conf.pm) (see [documentation](https://github.com/Ensembl/ensembl-production-imported/blob/main/docs/EGGeneric.md)).

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
| BRC4_genome_loader | creates an [Ensembl core database](http://www.ensembl.org/info/docs/api/core/index.html) from a set of flat files or adds ad-hoc (ie organellas) sequences to the existing core  | [BRC4_genome_loader](docs/BRC4_genome_loader.md) | [Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_loader_conf](lib/perl/Bio/EnsEMBL/Pipeline/PipeConfig/BRC4_genome_loader_conf.pm)
| BRC4_genome_dumper | | | | |
| | | | | |
| BRC4_genome_prepare | | | | |
| BRC4_addition_prepare | | | | |
| BRC4_genome_compare | | | | |
| | | | | |
| LoadGFF3 | | | | |
| LoadGFF3Batch | | | | |


### Scripts


## Various docs
See [docs](docs)

## TODO
Tests, tests, tests...

## Acknowledgements
Some of this code and documentation is inherited from the [EnsemblGenomes](https://github.com/EnsemblGenomes) and other [Ensembl](https://github.com/Ensembl) projects.
We appreciate the effort and time spent by developers of the [EnsemblGenomes](https://github.com/EnsemblGenomes) and [Ensembl](https://github.com/Ensembl) projects. 

Thank you!






