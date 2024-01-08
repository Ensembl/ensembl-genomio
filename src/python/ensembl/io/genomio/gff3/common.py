# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""A class that contains metadata to process GFF3 files:
- biotypes (supported, ignored)
"""

__all__ = ["GFFMeta"]

import logging
from typing import Any, Dict, List
from importlib_resources import files

from ensembl.io.genomio.utils.json_utils import get_json
import ensembl.io.genomio.data.gff3


class GFFMeta:
    """Class to share the list of feature types supported or ignored by the parser."""

    _biotypes: Dict[str, Dict[str, List[str]]] = {}

    def __init__(self) -> None:
        biotype_json = files(ensembl.io.genomio.data.gff3) / "biotypes.json"
        self._biotypes = get_json(biotype_json)
        logging.debug(f"Imported data from {biotype_json}")

    def get_biotypes(self, gene_type: str, supported: bool = True) -> List[str]:
        """Returns a list of biotypes supported or ignored.

        Args:
            gene_type: Gene type among "gene", "non_gene" or "transcript".
            supported: Fetch supported biotypes (otherwise returns the list of ignored biotypes).

        """
        if supported:
            return self._biotypes[gene_type]["supported"]
        return self._biotypes[gene_type]["ignored"]
