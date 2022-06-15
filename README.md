# ensembl-genomio
Pipelines to turn basic genomic data into Ensembl cores and back

A set of ehive pipelines:

* Genome dumper: dumps flat files from an Ensembl core database.
* Genome loader: creates an Ensembl core database from a set of flat files.

The files are in various formats (fasta, gff3, json) following BRC4 specifications.

## Get repo and install

```
git clone git@github.com:Ensembl/ensembl-genomio.git 
pip install ./ensembl-genomio

# test
python -c 'import ensembl.brc4.runnable.read_json'
```

Update your perl envs (if you need to)
```
export PERL5LIB=$(pwd)/ensembl-genomio/lib/perl:$PERL5LIB
export PATH=$(pwd)/ensembl-genomio/scripts:$PATH
```

### Optional installation
If you need to install "editable" python package use '-e' option
```
pip install -e ./ensembl-genomio
```

To install additional dependencies (e.g. `[doc]` or `[dev]`) provide `[<tag>]` string. I.e.
```
pip install -e ./ensembl-genomio[dev]
```

For the list of tags see `[project.optional-dependencies]` in [pyproject.toml](./pyproject.toml). 

### Docs generation
Install python part with the `[doc]` or `[dev]` tag.
Change into repo dir
Run doc build script.

```
git clone git@github.com:Ensembl/ensembl-genomio.git 
pip install -e ./ensembl-genomio[doc]

cd ./ensembl-genomio

# build docs
./scripts/setup/docs/build_sphinx_docs.sh
```



