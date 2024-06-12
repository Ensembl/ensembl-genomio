#!/usr/bin/sh

# Get the latest version of NCBI's datasets-cli software
DS_SOFTWARE_URL="https://api.github.com/repos/ncbi/datasets/releases/latest"
DATASETS_LATEST=`curl -s $DS_SOFTWARE_URL | grep browser_download_url | cut -d \" -f4 | grep linux-amd64.cli.package.zip`
DATASETS_RELEASE=`echo "$DATASETS_LATEST" | cut -d "/" -f 8`

# Export version info to bash file
echo "DATASETS_VERSION=$DATASETS_RELEASE"