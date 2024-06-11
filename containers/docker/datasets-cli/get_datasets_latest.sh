#!/usr/bin/bash

# Get the latest version of NCBI's datasets-cli software
DS_SOFTWARE_URL="https://api.github.com/repos/ncbi/datasets/releases/latest"

DATASETS_LATEST=`/usr/local/bin/curl -s $DS_SOFTWARE_URL | /usr/local/bin/grep browser_download_url | cut -d \" -f4 | /usr/local/bin/grep linux-amd64.cli.package.zip`

DATASETS_RELEASE=`echo "$DATASETS_LATEST" | cut -d "/" -f 8`
export DATASETS_RELEASE

# Get the latest 'jq' software
JQ_SOFTWARE_URL="https://api.github.com/repos/jqlang/jq/releases/latest"
JQ_LATEST=`/usr/local/bin/curl -s $JQ_SOFTWARE_URL | /usr/local/bin/grep browser_download_url | cut -d \" -f4 | /usr/local/bin/grep jq-linux-amd64`

cd /usr/local/bin/
wget -q $DATASETS_LATEST
unzip /usr/local/bin/linux-amd64.cli.package.zip
rm /usr/local/bin/linux-amd64.cli.package.zip
chmod +x /usr/local/bin/datasets /usr/local/bin/dataformat
wget $JQ_LATEST -O /usr/local/bin/jq
chmod +x /usr/local/bin/jq
