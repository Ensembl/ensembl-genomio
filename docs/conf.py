# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from datetime import datetime, UTC
import os
from pathlib import Path
import sys

import ensembl.utils

sys.path.insert(0, str(Path(__file__).parent / "_ext"))
sys.path.insert(0, str(Path(__file__).parent.parent / "src" / "python"))

# Project information
project = "ensembl-genomio"
author = "EMBL-European Bioinformatics Institute"
copyright = f"2016-{datetime.now(tz=UTC).year}, EMBL-European Bioinformatics Institute"  # noqa: A001

# General configuration
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "entrypoints_table",
]
language = "en"

# MyST settings
myst_enable_extensions = ["colon_fence", "substitution"]
myst_heading_anchors = 2

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
    "private-members": False,
    "undoc-members": False,
    "special-members": "__repr__",
}
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented"
suppress_warnings = ["autodoc.duplicate_object", "sphinx_autodoc_typehints.forward_reference"]
typehints_defaults = "comma"
typehints_document_rtype_none = False

# Napolean settings
napoleon_use_ivar = True

# Coverage settings
coverage_write_headline = False

# HTML output settings
html_theme = "pydata_sphinx_theme"
html_logo = "_static/ensembl_mark_white.png"
html_favicon = "_static/ensembl_favicon.png"
html_sourcelink_suffix = ""
html_last_updated_fmt = ""

html_title = "ensembl-genomio"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_extra_path = ["_static/htmlcov"]
html_js_files = [
    ("custom-icons.js", {"defer": "defer"}),
]

# Define the version and json_url for the version switcher
release = ensembl.utils.__version__
version_match = f"v{release}"
if os.environ.get("READTHEDOCS") or os.environ.get("CI"):
    json_url = "https://ensembl.github.io/ensembl-genomio/switcher.json"
else:
    json_url = "_static/switcher.json"

html_theme_options = {
    "header_links_before_dropdown": 4,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/Ensembl/ensembl-genomio",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/ensembl-genomio/",
            "icon": "fa-custom fa-pypi",
        },
    ],
    "logo": {
        "text": project,
    },
    "use_edit_page_button": False,
    "show_toc_level": 2,
    "navbar_align": "left",
    "show_version_warning_banner": True,
    "navbar_center": ["version-switcher", "navbar-nav"],
    "footer_start": ["copyright"],
    "footer_center": ["sphinx-version"],
    "secondary_sidebar_items": {
        "**/*": ["page-toc"],
        "coverage_report": [],
    },
    "switcher": {
        "json_url": json_url,
        "version_match": version_match,
    },
    "search_as_you_type": True,
    "navigation_with_keys": True,
}

html_sidebars = {
    "coverage_report": [],
}

html_context = {
    "github_user": "Ensembl",
    "github_repo": "ensembl-genomio",
    "github_version": "main",
    "doc_path": "docs",
}
