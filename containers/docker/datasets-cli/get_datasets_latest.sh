#!/usr/bin/sh

# Get the latest version of NCBI's datasets-cli software
DATASETS_LATEST="$1"

JQ_SOFTWARE_URL="https://api.github.com/repos/jqlang/jq/releases/latest"
JQ_LATEST=$(wget "$JQ_SOFTWARE_URL" -O - | grep browser_download_url | grep jq-linux-amd64 |  cut -f 4 -d '"')

# getting proper datasets binary
wget -q "$DATASETS_LATEST"
unzip linux-amd64.cli.package.zip
rm linux-amd64.cli.package.zip
chmod +x datasets dataformat

# getting jq util binary
wget "$JQ_LATEST" -O jq
chmod +x jq
