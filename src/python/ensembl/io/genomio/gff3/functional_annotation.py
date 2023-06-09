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
"""Simple representation of gene features functional annotation extracted from a GFF3 file."""

from os import PathLike
from pathlib import Path
import re
from typing import Any, Dict, List

from Bio.SeqFeature import SeqFeature

from ensembl.io.genomio.utils.json_utils import print_json


Annotation = Dict[str, Any]


class FunctionalAnnotations:
    """List of annotations extracted from a GFF3 file."""

    def __init__(self) -> None:
        self.annotations: List[Annotation] = []

    def add_feature(self, feature: SeqFeature, feat_type: str) -> None:
        """Append a feature object following the specifications.

        Args:
            feature: The SeqFeature to add to the list.
            feat_type: Feature type of the feature to store (e.g. gene, transcript, translation).

        """

        feature_object: Annotation = {"object_type": feat_type, "id": feature.id}

        # Description?
        if "product" in feature.qualifiers:
            description = feature.qualifiers["product"][0]
            if self._product_is_valid(description):
                feature_object["description"] = description

        if "Name" in feature.qualifiers and "description" not in feature_object:
            name = feature.qualifiers["Name"][0]

            # Exclude Name if it just a variant of the feature ID
            if feature.id not in name:
                feature_object["description"] = name

        # Synonyms?
        if "Name" in feature.qualifiers:
            feat_name = feature.qualifiers["Name"][0]
            if feat_name != feature.id:
                feature_object["synonyms"] = {"synonym": feat_name, "default": True}

        # is_pseudogene?
        if feature.type.startswith("pseudogen"):
            feature_object["is_pseudogene"] = True

        self.annotations.append(feature_object)

    @staticmethod
    def _product_is_valid(product: str) -> bool:
        """Returns True if the product name is valid, False otherwise.

        Args:
            product: A product name.

        """
        excluded_names = re.compile(
            r"^(uncharacterized|putative|hypothetical|predicted)"
            r"( uncharacterized)?"
            r" protein"
            r"( of unknown function)?"
            r"( \(fragment\))?$"
        )

        product_is_valid = not excluded_names.match(product.lower())
        return product_is_valid

    def _cleanup(self) -> None:
        """Returns the functional annotations list without putative product descriptions."""
        for feat in self.annotations:
            if "description" in feat and not self._product_is_valid(feat["description"]):
                del feat["description"]

    def to_json(self, out_path: PathLike) -> None:
        """Print out the current annotation list in a json file.

        Args:
            out_path: JSON file path where to write the data.

        """
        self._cleanup()
        print_json(Path(out_path), self.annotations)
