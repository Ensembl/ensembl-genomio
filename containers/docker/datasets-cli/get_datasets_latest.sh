#!/usr/bin/sh

# Get the latest version of NCBI's datasets-cli software
DS_SOFTWARE_URL="https://api.github.com/repos/ncbi/datasets/releases/latest"
DATASETS_LATEST=`curl -s $DS_SOFTWARE_URL | grep browser_download_url | cut -d \" -f4 | grep linux-amd64.cli.package.zip`
DATASETS_RELEASE=`echo "$DATASETS_LATEST" | cut -d "/" -f 8`

JQ_SOFTWARE_URL="https://api.github.com/repos/jqlang/jq/releases/latest"
JQ_LATEST=`curl -s $JQ_SOFTWARE_URL | grep browser_download_url | cut -d \" -f4 | grep jq-linux-amd64`

wget -q $DATASETS_LATEST
unzip linux-amd64.cli.package.zip
rm linux-amd64.cli.package.zip
chmod +x datasets dataformat
wget $JQ_LATEST -O jq
chmod +x jq
