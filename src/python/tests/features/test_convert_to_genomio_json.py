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
"""Unit testing of the ``convert_to_genomio_json`` package facade."""

from ensembl.io.genomio.features import convert_to_genomio_json
from ensembl.io.genomio.features.convert_to_genomio_json import args as convert_args
from ensembl.io.genomio.features.convert_to_genomio_json import (
    converters,
    document,
    registry,
    repeatmasker,
    trf,
)


def test_package_modules_expose_expected_entry_points() -> None:
    """Test ``convert_to_genomio_json`` package modules expose their owned entry points."""
    assert convert_args.parse_args is convert_to_genomio_json.parse_args
    assert converters.ConverterOptions is convert_to_genomio_json.ConverterOptions
    assert converters.FeatureConverter is convert_to_genomio_json.FeatureConverter
    assert registry.CONVERTERS_BY_LOGIC_NAME is convert_to_genomio_json.CONVERTERS_BY_LOGIC_NAME
    assert document.create_genomio_json is convert_to_genomio_json.create_genomio_json
    assert repeatmasker.parse_repeatmasker_output is convert_to_genomio_json.parse_repeatmasker_output
    assert trf.parse_trf_output is convert_to_genomio_json.parse_trf_output
