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
"""GFF parsing exceptions."""

__all__ = [
    "GeneSegmentError",
    "GFFParserError",
    "IgnoredFeatureError",
    "UnsupportedFeatureError",
]


class GFFParserError(Exception):
    """Error when parsing a GFF3 file."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.message = message


class GeneSegmentError(GFFParserError):
    """GFF3 gene segment parsing error."""


class IgnoredFeatureError(GFFParserError):
    """GFF3 feature can be ignored."""


class UnsupportedFeatureError(GFFParserError):
    """GFF3 feature is not supported."""
