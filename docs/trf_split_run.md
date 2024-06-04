# [trf_split_run.bash](https://github.com/Ensembl/ensembl-genomio/blob/main/scripts/trf_split_run.bash)

A trf wrapper with chunking support to be used with
[ensembl-production-imported DNAFeatures pipeline](https://github.com/Ensembl/ensembl-production-imported/tree/main/src/perl/Bio/EnsEMBL/EGPipeline/PipeConfig/DNAFeatures_conf.pm)
Compatible compatible with input/output format of trf invocation from [Bio::EnsEMBL::Analysis::Runnable::TRF](https://github.com/Ensembl/ensembl-analysis/blob/main/modules/Bio/EnsEMBL/Analysis/Runnable/TRF.pm).
And can be used as a hack to allow TRF stage to be accomplished at the cost of splitting
 long repeat into several adjacent ones  (with possible losses).

## Prerequisites
You should have [Biopython](https://biopython.org) installed and available in your environment.
You may check this with
```
python -c 'from Bio import SeqIO' || echo "no biopython" >> /dev/stderr
```

## Options
Use environment variable to control script run
* `DNA_FEATURES_TRF_SPLIT_NO_SPLITTING` -- set to `YES` to skip splitting stage
* `DNA_FEATURES_TRF_SPLIT_NO_TRF` -- set to `YES` to skip trf stage
* `DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE` -- chunk size [`1_000_000`]
* `DNA_FEATURES_TRF_SPLIT_SPLITTER_OPTIONS` -- for a finer control [`--n_seq 1 --chunk_tolerance 10 --chunk_size ${DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE}`] 
* `DNA_FEATURES_TRF_SPLIT_TRF_EXE` -- trf executable (or abs path to be used) [`trf`]
* `DNA_FEATURES_TRF_SPLIT_TRF_OPTIONS` -- additional options for TRF (like `-l 10`) []

## Usage examples
### A standalone run
```
trf_split_run.bash /writable/path_to/dna.fasta 2 5 7 80 10 40 500 -d -h
```

### As a substitution for the "trf_exe" to be used with the `TRF` stage of the `DNAFeatures` pipeline:
```
# change TRF "program" parameter, you should have ENSEMBL_ROOT_DIR env defined
tweak_pipeline.pl -url "$DNA_FEATURES_EHIVE_DB_URL" -tweak 'analysis[TRF].param[parameters_hash]={program=>"'${ENSEMBL_ROOT_DIR}'/ensembl-genomio/scripts/trf_split_run.bash"}'
```
```
# set environment variables if you need to, i.e.
export DNA_FEATURES_TRF_SPLIT_TRF_EXE=trf.4.09.1 
export DNA_FEATURES_TRF_SPLIT_TRF_OPTIONS='-l 10' # N.B. "l" correlated with the chunk size (-l chunk_size / 10^6)
export DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE=10_000_000'
```
```
runWorker.pl ... 
# or
beekeeper.pl ... -loop
```
```
# revert back (if you need to)
tweak_pipeline.pl -url $DNA_FEATURES_EHIVE_DB_URL  -tweak 'analysis[TRF].param[parameters_hash]={}'
```
```
# revert env settings
unset DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE
unset DNA_FEATURES_TRF_SPLIT_SPLITTER_OPTIONS
unset DNA_FEATURES_TRF_SPLIT_TRF_EXE
unset DNA_FEATURES_TRF_SPLIT_TRF_OPTIONS
#   script stages
unset DNA_FEATURES_TRF_SPLIT_NO_SPLITTING
unset DNA_FEATURES_TRF_SPLIT_NO_TRF
```

