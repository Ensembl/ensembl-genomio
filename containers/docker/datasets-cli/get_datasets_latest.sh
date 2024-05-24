#!/usr/bin/sh

# Get the latest version of NCBI's datasets-cli software
DS_SOFTWARE_URL="https://api.github.com/repos/ncbi/datasets/releases/latest"
export DATASETS_LATEST=`curl -s $DS_SOFTWARE_URL | grep browser_download_url | cut -d \" -f4 | egrep linux-amd64.cli.package.zip`
export DATASETS_RELEASE=`echo "$DATASETS_LATEST" | cut -d "/" -f 8`

JQ_SOFTWARE_URL="https://api.github.com/repos/jqlang/jq/releases/latest"
export JQ_LATEST=`curl -s $JQ_SOFTWARE_URL | grep browser_download_url | cut -d \" -f4 | egrep jq-linux-amd64`
