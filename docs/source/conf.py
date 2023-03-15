# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import datetime

from mock import Mock as MagicMock


class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()


MOCK_MODULES = []
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

sys.path.insert(0, os.path.abspath("../../src/python"))

print(sys.executable)

# -- Project information -----------------------------------------------------

project = "ensembl-genomio"
author = "Ensembl Metazoa"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = "0.1"
# The full version, including alpha/beta/rc tags.
release = "0.1"

copyright_owner = "EMBL-European Bioinformatics Institute"
copyright_dates = "[2016-%d]" % datetime.datetime.now().year
copyright = copyright_dates + " " + copyright_owner

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

source_suffix = {".rst": "restructuredtext"}

# The master toctree document.
master_doc = "index"

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
]

# A list of ignored prefixes for module index sorting.
modindex_common_prefix = ["ensembl.brc4.runnable"]

# Defining autodoc functionality
autodoc_default_options = {
    "member-order": "alphabetical",
    "undoc-members": False,
    "exclude-members": "__weakref__",
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = "agogo"
html_theme_options = {
    "bodyfont": "Garamond, Arial, serif",
    "headerfont": "Arial, Helvetica, serif",
    "headerlinkcolor": "#33d6ff",
    "pagewidth": "70em",
    "documentwidth": "50em",
    "rightsidebar": True,
    "bgcolor": "#009999",
    "headerbg": "#009999",
    "footerbg": "#e6fff9",
    "linkcolor": "green",
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "ensembl-genomio", "Ensembl Genomio Base Library Documentation", [author], 1)]

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "ensembl-genomio.tex", "Ensembl-genomio Base Library Documentation", [author], "manual"),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "ensembl-genomio", "Ensembl-genomio Base Library Documentation", [author], 1)]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "ensembl-genomio",
        "Ensembl-genomio Library Documentation",
        author,
        "ensembl-genomio",
        "Ensembl-genomio Base Library.",
        "Miscellaneous",
    ),
]
