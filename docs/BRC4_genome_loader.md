## BRC4_genome_loader
### Module: [Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_loader_conf](https://github.com/Ensembl/ensembl-genomio/blob/main/src/perl/Bio/EnsEMBL/Pipeline/PipeConfig/BRC4_genome_loader_conf.pm)

Creates an [Ensembl core database](http://www.ensembl.org/info/docs/api/core/index.html)
from a set of flat files
or adds ad-hoc (i.e. organellas) sequences to the existing core.

### Prerequisites
A registry file with the locations of the core database server(s) and the production database (or `-production_db $PROD_DB_URL` specified).

### How to run
```
# REP_LIB_OPT= # or whatever

init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::NOTDNAFeatures_conf \
  $($CMD details script) \
  -registry $REG_FILE \
  -production_db "$($PROD_SERVER details url)""$PROD_DBNAME" \
  -hive_force_init 1\
  -pipeline_tag "_${SPECIES_TAG}" \
  -pipeline_dir $OUT_DIR \
  -species $SPECIES \
  -redatrepeatmasker 0 \
  -always_use_repbase 1 \
  -repeatmasker_timer '10H' \
  $REP_LIB_OPT \
  -repeatmasker_repbase_species "$REPBASE_SPECIES_NAME" \
  -max_seq_length 300000 \
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


### Parameters / Options

| option | default value |  meaning |
| - | - | - |
| `-division` | | division (intersection with registry) to be processed
| `-species` |  | species to process, several `-species ` options are possible
| `-pipeline_dir` | | directory to store results to



#### Notes
2nd mode for adding



### Input data


### Parts


#### Runnables


#### Configuration files
xref map


