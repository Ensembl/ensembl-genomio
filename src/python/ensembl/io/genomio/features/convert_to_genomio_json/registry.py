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
"""Explicit registry of supported feature output converters."""

from ensembl.io.genomio.features.convert_to_genomio_json.repeatmasker import (
    RepeatMaskerConverter,
    RepeatMaskerCustomConverter,
    RepeatMaskerRepbaseConverter,
)
from ensembl.io.genomio.features.convert_to_genomio_json.trf import TrfConverter

__all__ = [
    "CONVERTERS_BY_LOGIC_NAME",
    "TOP_LEVEL_CONVERTERS",
]

CONVERTERS_BY_LOGIC_NAME = {
    RepeatMaskerCustomConverter.analysis_logic_name: RepeatMaskerCustomConverter,
    RepeatMaskerRepbaseConverter.analysis_logic_name: RepeatMaskerRepbaseConverter,
    TrfConverter.analysis_logic_name: TrfConverter,
}

TOP_LEVEL_CONVERTERS = (
    TrfConverter,
    RepeatMaskerConverter,
)
