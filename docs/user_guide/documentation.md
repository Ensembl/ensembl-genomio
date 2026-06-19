# Documentation

This library's documentation is generated using _sphinx_'s PyData theme. It is automatically handled by a GitHub Actions workflow, but if you want to generate a local instance for development or testing, you can create a virtual environment, clone the repository, and install GenomIO with the `docs` tag to install the sphinx framework dependencies as shown below:

```bash
python -m venv <VIRTUAL_ENVIRONMENT_NAME>
source <VIRTUAL_ENVIRONMENT_NAME>/bin/activate
git clone https://github.com/Ensembl/ensembl-genomio.git
cd ensembl-genomio
pip install -e .[docs]
```

Once the installation is completed, generating the documentation is very easy:

```bash
make docs
```

This will generate a folder inside `docs` called `_build`. You can open the `docs/_build/html/index.html` file with your preferred browser to explore your newly created instance of the documentation.

If you want to include coverage information, you will need to generate the coverage information first:

```bash
make coverage
make docs
```

> :warning: Bear in mind that in order to generate the coverage report, all the unit tests will need to be run first, so additional compute resources may be required. You will also need to add the `cicd` tag to the installation step described above:
>
> ```bash
> pip install -e .[docs,cicd]
> ```

For full information on the framework used to generate the documentation, visit [sphinx-doc.org](https://www.sphinx-doc.org/en/master/index.html) and [pydata-sphinx-theme](https://pydata-sphinx-theme.readthedocs.io/).
