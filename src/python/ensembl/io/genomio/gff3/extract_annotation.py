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

__all__ = [
    "Annotation",
    "DuplicateIdError",
    "MissingParentError",
    "AnnotationError",
    "FunctionalAnnotations",
]

import logging
from os import PathLike
from pathlib import Path
import re
from typing import Any, Dict, List, Optional

from Bio.SeqFeature import SeqFeature

from ensembl.io.genomio.utils.json_utils import print_json


Annotation = Dict[str, Any]

_PARENTS = {
    "transcript": "gene",
    "translation": "transcript",
}


class DuplicateIdError(Exception):
    """Trying to add a feature with an ID already in use."""


class MissingParentError(Exception):
    """Trying to add a feature without an expected parent."""


class AnnotationError(Exception):
    """If anything wrong happens when recording annotations."""


class FunctionalAnnotations:
    """List of annotations extracted from a GFF3 file."""

    def __init__(self) -> None:
        self.annotations: List[Annotation] = []

        # Annotated features
        # Under each feature, each dict's key is a feature ID
        self.features: Dict[str, Dict[str, Annotation]] = {
            "gene": {},
            "transcript": {},
            "translation": {},
            "transposable_element": {},
        }
        # Keep parent info: key is the feature ID, value is the parent ID
        self.parents: Dict[str, Dict[str, str]] = {
            "gene": {},
            "transcript": {},
        }

    def get_features(self, feat_type: str) -> Dict[str, Annotation]:
        """Get all feature annotations for the requested type."""
        try:
            return self.features[feat_type]
        except KeyError as err:
            raise KeyError(f"No such feature type {feat_type}") from err

    def add_parent_link(self, parent_type: str, parent_id: str, child_id: str) -> None:
        """Record a parent-child IDs relationship for a given parent biotype."""
        features = self.get_features(parent_type)
        if parent_id not in features:
            raise MissingParentError(f"Parent {parent_type}:{parent_id} not found for {child_id}")
        self.parents[parent_type][child_id] = parent_id

    def get_parent(self, parent_type: str, child_id: str) -> str:
        """Returns the parent ID of a given child for a given parent biotype."""
        try:
            parents = self.parents[parent_type]
        except KeyError as err:
            raise KeyError(f"Unsupported parent type {parent_type}") from err

        parent_id = parents.get(child_id)
        if parent_id is None:
            raise MissingParentError(f"Can't find {parent_type} parent for {child_id}")
        return parent_id

    def add_feature(
        self,
        feature: SeqFeature,
        feat_type: str,
        parent_id: Optional[str] = None,
        all_parent_ids: Optional[List[str]] = None,
    ) -> None:
        """Add annotation for a feature of a given type. If a parent_id is provided, record the relatioship.

        Args:
            feature: The feature to create an annotation.
            feat_type: Type of the feature to annotate.
        """
        if all_parent_ids is None:
            all_parent_ids = []
        features = self.get_features(feat_type)
        if feature.id in features:
            raise AnnotationError(f"Feature {feat_type} ID {feature.id} already added")

        feature_object = self._generic_feature(feature, feat_type, all_parent_ids)
        self.features[feat_type][feature.id] = feature_object

        if parent_id:
            if feat_type in _PARENTS:
                parent_type = _PARENTS[feat_type]
                self.add_parent_link(parent_type, parent_id, feature.id)
            else:
                raise AnnotationError(f"No parent possible for {feat_type} {feature.id}")

    def _generic_feature(
        self, feature: SeqFeature, feat_type: str, parent_ids: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """Create a feature object following the specifications.

        Args:
            feature: The SeqFeature to add to the list.
            feat_type: Feature type of the feature to store (e.g. gene, transcript, translation).

        """
        if parent_ids is None:
            parent_ids = []

        feature_object: Annotation = {"object_type": feat_type, "id": feature.id}

        # Description?
        for qname in ("description", "product"):
            if qname in feature.qualifiers:
                description = feature.qualifiers[qname][0]
                all_ids = list(parent_ids)
                all_ids.append(feature.id)
                if self.product_is_informative(description, feat_ids=all_ids):
                    feature_object["description"] = description
                    break
                logging.debug(f"Non informative description for {feature.id}: {description}")

        # Synonyms?
        if "Name" in feature.qualifiers:
            feat_name = feature.qualifiers["Name"][0]
            if feat_name != feature.id:
                feature_object["synonyms"] = [{"synonym": feat_name, "default": True}]

        # is_pseudogene?
        if feature.type.startswith("pseudogen"):
            feature_object["is_pseudogene"] = True

        return feature_object

    @staticmethod
    def product_is_informative(product: str, feat_ids: Optional[List[str]] = None) -> bool:
        """Returns True if the product name contains informative words, False otherwise.

        It is considered uninformative when the description contains words such as "hypothetical" or
        or "putative". If a feature IDs are provided, consider it uninformative as well (we do not want
        descriptions to be just the ID).

        Args:
            product: A product name.
            feat_ids: List of feature ID.

        """
        if feat_ids is None:
            feat_ids = []
        non_informative_words = [
            "hypothetical",
            "putative",
            "uncharacterized",
            "unspecified",
            "unknown",
            r"(of )?unknown function",
            "conserved",
            "predicted",
            "fragment",
            "product",
            "function",
            "protein",
            "transcript",
            "gene",
            "RNA",
            r"(variant|isoform)( X?\d+)?",
        ]
        non_informative_re = re.compile(r"|".join(non_informative_words), re.IGNORECASE)

        # Remove all IDs they are in the description
        if feat_ids:
            logging.debug(f"Filter out {feat_ids} from {product}")
            try:
                for feat_id in feat_ids:
                    feat_id_re = re.compile(feat_id, re.IGNORECASE)
                    product = re.sub(feat_id_re, "", product)
            except TypeError as err:
                raise TypeError(f"Failed to search {feat_id_re} in '{product}'") from err

        # Remove punctuations
        punct_re = re.compile(r"[,;: _()-]+")
        product = re.sub(punct_re, " ", product)

        # Then remove non informative words
        product = re.sub(non_informative_re, " ", product)

        # Anything (informative) left?
        empty_re = re.compile(r"^[ ]*$")
        return not bool(empty_re.match(product))

    def transfer_descriptions(self) -> None:
        """Transfers the feature descriptions in 2 steps:
        - from translations to transcripts (if the transcript description is empty)
        - from transcripts to genes (same case)

        """
        self._transfer_description_up("translation")
        self._transfer_description_up("transcript")

    def _transfer_description_up(self, child_feature: str) -> None:
        """Transfer descriptions from all feature of a given type, up to their parent.

        Args:
            child_feature: Either "translation" (transfer to transcript) or "transcript" (to gene).

        """
        children_features = self.get_features(child_feature)
        parent_type = _PARENTS[child_feature]
        parent_features = self.get_features(parent_type)

        # Transfer description from children to their parent
        for child_id, child in children_features.items():
            child_description = child.get("description")
            if child_description is not None:
                # Check parent
                parent_id = self.get_parent(parent_type, child_id)
                parent = parent_features[parent_id]
                parent_description = parent.get("description")
                if parent_description is None:
                    parent["description"] = child_description

    def _to_list(self):
        all_list = []
        for feat_dict in self.features.values():
            all_list += feat_dict.values()
        return all_list

    def to_json(self, out_path: PathLike) -> None:
        """Print out the current annotation list in a json file.

        Args:
            out_path: JSON file path where to write the data.

        """
        self.transfer_descriptions()
        feats_list = self._to_list()
        print_json(Path(out_path), feats_list)
