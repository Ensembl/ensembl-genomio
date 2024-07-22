# Installation

GenomIO is now publicly available in [PyPI](https://pypi.org), so you can easily install it with your preferred Python package manager, e.g.

```bash
pip install ensembl-genomio
```

## API setup and installation

To run also the pipelines present in the repository, you will need to have Perl 5.26 and Bioperl 1.6.9+ installed in your system, and then clone the GitHub repository together with some other repositories:

```
git clone https://github.com/Ensembl/ensembl-genomio
git clone https://github.com/Ensembl/ensembl-analysis -b dev/hive_master
git clone https://github.com/Ensembl/ensembl-production
git clone https://github.com/Ensembl/ensembl-hive
git clone https://github.com/Ensembl/ensembl-taxonomy
git clone https://github.com/Ensembl/ensembl-orm
```

### Documentation

Documentation for Ensembl-genomio generated using _mkdocs_. For full information visit [mkdocs.org](https://www.mkdocs.org).

Useful commands:
* `mkdocs new [dir-name]` - Create a new project
* `mkdocs serve` - Start the live-reloading docs server
* `mkdocs build` - Build the documentation site
* `mkdocs -h` - Print help message and exit
