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
"""Sphinx configuration for the ensembl-genomio documentation.

All shared defaults live in :mod:`ensembl.utils.docs`; this file only sets project-specific values and any
local overrides.
"""

from pathlib import Path

import ensembl.io.genomio
from ensembl.utils.docs import configure

coverage_root = Path(__file__).parent / "reports"
configure(
    globals(),
    project="ensembl-genomio",
    repo_url="https://github.com/Ensembl/ensembl-genomio",
    release=ensembl.io.genomio.__version__,
    docs_base_url="https://ensembl.github.io/ensembl-genomio",
    coverage_root=coverage_root if coverage_root.exists() else None,
    include_entrypoints=True,
    add_pypi_icon=True,
)
