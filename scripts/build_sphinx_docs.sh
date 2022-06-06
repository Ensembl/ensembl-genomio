#!/bin/bash

sphinx-apidoc -Mf --implicit-namespaces -o docs/source/ lib/python/ensembl
cd docs
make clean && make html
