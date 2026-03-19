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
"""Unit testing of `ensembl.io.genomio.utils.json_utils` module."""

from email.mime import text
import json
from pathlib import Path

import pytest

from ensembl.io.genomio.utils import json_utils


def test_get_json_loads_data(tmp_path: Path):
    p = tmp_path / "in.json"
    payload = {"b": 2, "a": 1, "nested": {"x": [1, 2, 3]}}
    p.write_text(json.dumps(payload), encoding="utf-8")

    out = json_utils.get_json(p)
    assert out == payload


def test_get_json_passes_parse_int_kwarg(tmp_path):
    p = tmp_path / "data.json"
    p.write_text('{"value": 123}', encoding="utf-8")

    def parse_int_as_str(x: str) -> str:
        return f"INT:{x}"

    out = json_utils.get_json(p, parse_int=parse_int_as_str)

    assert out["value"] == "INT:123"


def test_print_json_uses_defaults_sort_keys_and_indent(tmp_path: Path):
    p = tmp_path / "out.json"

    json_utils.print_json(p, {"b": 2, "a": 1})

    text = p.read_text(encoding="utf-8")

    # Default indent=4 should introduce newline + 4 spaces before keys
    assert "\n    " in text

    # Default sort_keys=True means "a" should appear before "b"
    assert text.index('"a"') < text.index('"b"')

    assert json.loads(text) == {"a": 1, "b": 2}


def test_print_json_respects_overrides_for_sort_keys_and_indent(tmp_path: Path):
    p = tmp_path / "out.json"

    json_utils.print_json(p, {"b": 2, "a": 1}, sort_keys=False, indent=2)

    text = p.read_text(encoding="utf-8")

    assert "\n  " in text
    assert "\n    " not in text
    assert text.index('"b"') < text.index('"a"')
    assert json.loads(text) == {"a": 1, "b": 2}


def test_print_json_allows_non_default_serialization_options(tmp_path: Path):
    p = tmp_path / "out.json"

    # separators removes spaces after commas/colons
    json_utils.print_json(p, {"a": 1, "b": 2}, separators=(",", ":"))

    text = p.read_text(encoding="utf-8")
    assert '": ' not in text
    assert json.loads(text) == {"a": 1, "b": 2}
