# Building the documentation locally

This library's documentation is generated using _sphinx_'s PyData theme. It is automatically handled by a GitHub Actions workflow, but if you want to generate a local instance for development or testing, you can create a virtual environment, clone the repository, and install this library with the `docs` tag to install the sphinx framework dependencies as shown below:

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
make coverage && make docs
```

```{warning}
Bear in mind that in order to generate the coverage report, all the unit tests will need to be run first, so additional compute resources may be required. You will also need to add the `cicd` tag to the installation step described above, i.e. `pip install -e .[cicd,docs]`.
```

For full information on the framework used to generate the documentation, visit [sphinx-doc.org](https://www.sphinx-doc.org/en/master/index.html) and [pydata-sphinx-theme](https://pydata-sphinx-theme.readthedocs.io/).

## Previewing locally

Open `docs/_build/html/index.html` directly in your browser for a quick look at the rendered content. However, some features require the docs to be served over HTTP rather than loaded from the filesystem, like the version switcher: browsers block this when the page is opened as a `file://` URL, so the version dropdown will appear empty. Serve the build locally to work around this:

```bash
python -m http.server 8000 --directory docs/_build/html
```

Then open `http://localhost:8000` in your browser. The switcher will populate and highlight the current version correctly (see [below](#generating-the-version-switcher-locally) for more information).

### Generating the version switcher locally

The `update_docs_switcher` CLI generates the `switcher.json` expected by the theme. For a local preview, first generate the documentation with `make docs` as indicated above, then run the following command with the version you want to test against:

```bash
update_docs_switcher --base-url http://localhost:8000 --version 1.0.0 docs/_build/html/_static/switcher.json
```

Finally, serve the build locally as shown above.

```{note}
The `release` value set in `docs/conf.py` must be an exact, case-sensitive match to one of the `version`
fields in `switcher.json` — including any `v` prefix. A mismatch causes the switcher to render without a
highlighted selection. Check the browser console for a warning from the theme if the dropdown appears empty.
```
