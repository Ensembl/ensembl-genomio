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

__all__ = [
    "GFFMeta"
]

import logging
from typing import Any, Dict, List
from importlib_resources import files

from ensembl.io.genomio.utils.json_utils import get_json
import ensembl.io.genomio.data.gff3


class GFFMeta:
    """Heritable class to share the list of feature types supported or ignored by the parser"""
    _biotypes = {}

    @classmethod
    def _load_biotypes(cls) -> None:
        biotype_json = files(ensembl.io.genomio.data.gff3) / "biotypes.json"
        logging.debug(f"Import data from {biotype_json}")
        cls._biotypes = get_json(biotype_json)

    @classmethod
    def get_biotypes(cls) -> Dict[str, Any]:
        if not cls._biotypes:
            cls._load_biotypes()
        return cls._biotypes

    @classmethod
    def gene_types(cls) -> List[str]:
        if not cls._biotypes:
            cls._load_biotypes()
        return cls._biotypes["gene"]["supported"]

    @classmethod
    def ignored_gene_types(cls) -> List[str]:
        if not cls._biotypes:
            cls._load_biotypes()
        return cls._biotypes["gene"]["ignored"]

    @classmethod
    def non_gene_types(cls) -> List[str]:
        if not cls._biotypes:
            cls._load_biotypes()
        return cls._biotypes["non_gene"]["supported"]

    @classmethod
    def transcript_types(cls) -> List[str]:
        if not cls._biotypes:
            cls._load_biotypes()
        return cls._biotypes["transcript"]["supported"]

    @classmethod
    def ignored_transcript_types(cls) -> List[str]:
        if not cls._biotypes:
            cls._load_biotypes()
        return cls._biotypes["transcript"]["ignored"]
