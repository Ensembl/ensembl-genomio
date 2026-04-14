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
"""GFF3 features."""

from __future__ import annotations

__all__ = [
    "GFFSeqFeature",
]

from Bio.SeqFeature import SeqFeature, Location


class GFFSeqFeature(SeqFeature):
    """Extends `Bio.SeqFeature.SeqFeature` with sub_features, to be used for typing."""

    def __init__(
        self,
        location: Location | None = None,
        *,
        type: str = "",  # pylint: disable=W0622
        id: str = "<unknown id>",  # pylint: disable=W0622
        qualifiers: dict | None = None,
        sub_features: list[GFFSeqFeature] | None = None,
    ):
        super().__init__(location, type=type, id=id, qualifiers=qualifiers)
        if sub_features is None:
            sub_features = []
        self.sub_features = sub_features

    @classmethod
    def cast(cls, feat: SeqFeature) -> GFFSeqFeature:
        """Cast a SeqFeature to a GFFSeqFeature."""
        feat.__class__ = cls
        if not hasattr(feat, "sub_features"):
            feat.sub_features = []  # type: ignore[attr-defined]
        return feat  # type: ignore[return-value]
