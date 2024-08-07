#!/usr/bin/env bash

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# a trf wrapper with chunking support to be used with DNAFeatures pipeline
#   compatible with  input/output format of trf invocation from  Bio::EnsEMBL::Analysis::Runnable::TRF
#       (https://github.com/Ensembl/ensembl-analysis/blob/main/modules/Bio/EnsEMBL/Analysis/Runnable/TRF.pm)
#   can be used as a hack to allow TRF stage to be accomplished at the cost of splitting long repeat into several adjacent ones
#   (with possible losses)
#   N.B. you should have Biopython installed and available in your environmnent. You may check this with:
#     python -c 'from Bio import SeqIO' || echo "no biopython" >> /dev/stderr
# usage:
#   use env vairable to control script run:
#    * "DNA_FEATURES_TRF_SPLIT_NO_SPLITTING" -- set to "YES" to skip splitting stage
#    * "DNA_FEATURES_TRF_SPLIT_NO_TRF" -- set to "YES" to skip trf stage
#    * "DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE" -- chunk size [1_000_000]
#    * "DNA_FEATURES_TRF_SPLIT_SPLITTER_OPTIONS" -- for a finer control [--n_seq 1 --chunk_tolerance 10 --chunk_size ${DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE}] 
#    * "DNA_FEATURES_TRF_SPLIT_TRF_EXE" -- trf executrable (or abs path to be used) [trf]
#    * "DNA_FEATURES_TRF_SPLIT_TRF_OPTIONS" -- addtitional options for TRF (like "-l 10", N.B. "-l chunk_size / 10^6") []
#
# examples
#  * standalone run
#    trf_split_run.bash /writable/path_to/dna.fasta 2 5 7 80 10 40 500 -d -h -l 10
#
#  * As a substitution for the "trf_exe" to be used with the `TRF` stage of the `DNAFeatures` pipeline:
#    # change TRF "program" parameter, you should have ENSEMBL_ROOT_DIR env defined
#    tweak_pipeline.pl -url "$DNA_FEATURES_EHIVE_DB_URL" -tweak 'analysis[TRF].param[parameters_hash]={program=>"'${ENSEMBL_ROOT_DIR}'/ensembl-genomio/scripts/trf_split_run.bash"}'
#    # set envitonment variables if you need to, i.e.
#    export DNA_FEATURES_TRF_SPLIT_TRF_EXE=trf.4.09.1 
#    runWorker.pl ... or beekeeper.pl ...
#
#    # revert back
#    tweak_pipeline.pl -url $DNA_FEATURES_EHIVE_DB_URL  -tweak 'analysis[TRF].param[parameters_hash]={}'
#
#    # env settings
#    unset DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE
#    unset DNA_FEATURES_TRF_SPLIT_SPLITTER_OPTIONS
#    unset DNA_FEATURES_TRF_SPLIT_TRF_EXE
#    unset DNA_FEATURES_TRF_SPLIT_TRF_OPTIONS
#
#    # script stages
#    unset DNA_FEATURES_TRF_SPLIT_NO_SPLITTING
#    unset DNA_FEATURES_TRF_SPLIT_NO_TRF


set -o errexit
set -o pipefail

# CONSTANTS
DEFAULT_CHUNK_SIZE=1_000_000
DEFAULT_CHUNK_TOLERANCE=10 # percent (int)
DEFAULT_CHUNK_SFX=TRFSPLTRN # do not use spaces or regexp special characters
DEFAULT_TRF_EXE=trf

# GETTING INPUT ARGUMENTS
ALL_ARGS=( "$@" )
FILE_NAME="${ALL_ARGS[@]:0:1}"
PARAMS="${ALL_ARGS[@]:1}"
PARAMS_FOR_NAME="${ALL_ARGS[@]:1:7}"

OUTPUT_FILE_NAME="$FILE_NAME"."$(echo $PARAMS_FOR_NAME | perl -pe 'chomp; s/\s+/./g')".dat
PARTS_DIR="$FILE_NAME"."$(echo $PARAMS_FOR_NAME | perl -pe 'chomp; s/\s+/./g')".dat.parts
echo "storing chnuks and partial results to $PARTS_DIR " >> /dev/stderr


# GETTING SCRIPT NAME
SCRIPT="$(which $0)"
echo "trf_split_run is called as $SCRIPT" >> /dev/stderr


# SLITTING
if [ -z "$DNA_FEATURES_TRF_SPLIT_NO_SPLITTING" -o x"$DNA_FEATURES_TRF_SPLIT_NO_SPLITTING" != x"YES" ]; then
  SPLITTER_SCRIPT="fasta_chunk"
  echo "using chunk_fasta sctipt: $SPLITTER_SCRIPT" >> /dev/stderr 

  SPLITTER_OPTIONS="--n_seq 1 --chunk_tolerance ${DEFAULT_CHUNK_TOLERANCE}"
  if [ -n "$DNA_FEATURES_TRF_SPLIT_SPLITTER_OPTIONS" ]; then
    echo "using DNA_FEATURES_TRF_SPLIT_SPLITTER_OPTIONS '$DNA_FEATURES_TRF_SPLIT_SPLITTER_OPTIONS' for splitting" >> /dev/stderr
    SPLITTER_OPTIONS="$DNA_FEATURES_TRF_SPLIT_SPLITTER_OPTIONS"
  else
    CHUNK_SIZE="${DEFAULT_CHUNK_SIZE}"
    if [ -n "$DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE" ]; then
      echo "using DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE '$DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE' for splitting" >> /dev/stderr
      CHUNK_SIZE="$DNA_FEATURES_TRF_SPLIT_SPLITTER_CHUNK_SIZE"
    fi
    SPLITTER_OPTIONS="$SPLITTER_OPTIONS --chunk_size $CHUNK_SIZE"
  fi
  echo "using spitter options: $SPLITTER_OPTIONS" >> /dev/stderr

  SPLIT_CMD="'$SPLITTER_SCRIPT' $SPLITTER_OPTIONS --chunk_sfx '${DEFAULT_CHUNK_SFX}' --add_offset --individual_out_dir '$PARTS_DIRi' --out part --fasta_dna '$FILE_NAME'"
  echo "Running \"$SPLIT_CMD\"" >> /dev/stderr
  "$SPLITTER_SCRIPT" $SPLITTER_OPTIONS --chunk_sfx "${DEFAULT_CHUNK_SFX}" --add_offset --individual_out_dir "$PARTS_DIR" --out part --fasta_dna "$FILE_NAME"
fi # DNA_FEATURES_TRF_SPLIT_NO_SPLITTING

# APPENDING OPTIONS
if [ -n "DNA_FEATURES_TRF_SPLIT_TRF_OPTIONS" ]; then
  echo "appending '$DNA_FEATURES_TRF_SPLIT_TRF_OPTIONS' to trf params" >> /dev/stderr
  PARAMS="$PARAMS $DNA_FEATURES_TRF_SPLIT_TRF_OPTIONS"
fi

# RUN TRF
if [ -z "$DNA_FEATURES_TRF_SPLIT_NO_TRF" -o x"$DNA_FEATURES_TRF_SPLIT_NO_TRF" != x"YES" ]; then

  if [ -z "$DNA_FEATURES_TRF_SPLIT_TRF_EXE" ]; then
    DNA_FEATURES_TRF_SPLIT_TRF_EXE="${DEFAULT_TRF_EXE}"
  fi
  echo "using DNA_FEATURES_TRF_SPLIT_TRF_EXE '$DNA_FEATURES_TRF_SPLIT_TRF_EXE' with '$PARAMS' for trf" >> /dev/stderr

  pushd "$PARTS_DIR"
    ls -1 part.*.fa |
      sort -k2,2n -k3,3n -t . |
      xargs -r -n 1 -I XXX sh -c "
        echo running '\"$DNA_FEATURES_TRF_SPLIT_TRF_EXE\" \"XXX\" $PARAMS' >> /dev/stderr
        \"$DNA_FEATURES_TRF_SPLIT_TRF_EXE\" \"XXX\" $PARAMS 2 >> /dev/stderr || [ $(($? % 256)) -eq 0 ]
      "
  popd

fi # DNA_FEATURES_TRF_SPLIT_NO_TRF

# MERGE
echo "doing merge step for '$PARTS_DIR'" >> /dev/stderr
pushd "$PARTS_DIR"
  ls -1 part.*.fa.*.dat |
    sort -k2,2n -k3,3n -t . |
    xargs -n 1 -I XXX cat XXX |
    perl -pe 's/^(Sequence:\s+[^\s]+)(?:_'"${DEFAULT_CHUNK_SFX}"'_[^\s]*_off_?)(\d+)\s.*\s*/$1 $2\n/' |
    awk '/^Sequence:/ {offset = 0; if (!seen[$2]) {print $1, $2}; seen[$2] = 1; offset = $3}
         /^[0-9]+/ {$1 += offset; $2 += offset; print}' |
    cat > "$OUTPUT_FILE_NAME"
popd

exit 0
