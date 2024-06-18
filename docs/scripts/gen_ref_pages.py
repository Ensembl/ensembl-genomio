# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Generate the code reference pages for mkdocs."""

from pathlib import Path

import mkdocs_gen_files


nav = mkdocs_gen_files.Nav()

root = Path("src/python").resolve()
ensembl_path = root / "ensembl"
for py_path in sorted(ensembl_path.rglob("*.py")):
    # Get the relative module path and corresponding documentation paths
    module_path = py_path.relative_to(root)
    doc_path = module_path.with_suffix(".md")
    full_doc_path = Path("reference", doc_path)
    # Get all the parts of the module path without the ".py" extension
    parts = tuple(module_path.with_suffix("").parts)
    # Drop "__init__" file from the path components as well (if present)
    if parts[-1] == "__init__":
        parts = parts[:-1]
        doc_path = doc_path.with_name("index.md")
        full_doc_path = full_doc_path.with_name("index.md")
    # Add markdown file path with its index tree
    nav[parts] = doc_path.as_posix()
    # Populate the markdown file with the doc stub of this module
    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        identifier = ".".join(parts)
        fd.write(f"::: {identifier}\n")
    # Correct the path
    mkdocs_gen_files.set_edit_path(full_doc_path, module_path)

with mkdocs_gen_files.open("reference/summary.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
