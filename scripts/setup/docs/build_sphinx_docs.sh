#!/bin/bash

cd docs
sphinx-apidoc -Mf --implicit-namespaces -o source ../lib/python/ensembl

SPHINX_MAKEFILE="../scripts/setup/docs/Makefile"
make -f $SPHINX_MAKEFILE clean && make -f $SPHINX_MAKEFILE html text
