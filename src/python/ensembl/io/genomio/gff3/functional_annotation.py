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
from typing import Any, Dict, List, Optional

from Bio.SeqFeature import SeqFeature

from ensembl.io.genomio.utils.json_utils import print_json


Annotation = Dict[str, Any]


class DuplicateIdError(Exception):
    """Trying to add a feature with an ID already in use."""


class MissingParentError(Exception):
    """Trying to add a feature without an expected parent."""


class FunctionalAnnotations:
    """List of annotations extracted from a GFF3 file."""

    def __init__(self) -> None:
        self.annotations: List[Annotation] = []

        # Annotated features
        self.genes: Dict[str, Annotation] = {}
        self.transcripts: Dict[str, Annotation] = {}
        self.translations: Dict[str, Annotation] = {}
        self.transposable_elements: Dict[str, Annotation] = {}
        # Keep parent info: key is the feature ID, value is the parent ID
        self.parent: Dict[str, str] = {}

    def add_parent(self, parent_type: str, parent_id: str, child_id: str) -> None:
        """Record a parent-child IDs relationship for a given parent biotype."""
        if parent_type in ("gene", "transcript"):
            self.parent[f"{parent_type}-{child_id}"] = parent_id
        else:
            raise MissingParentError(f"Unsupported parent type {parent_type}")

    def get_parent(self, parent_type: str, child_id: str) -> str:
        """Returns the parent ID of a given child for a given parent biotype."""
        if parent_type in ("gene", "transcript"):
            parent_id = self.parent.get(f"{parent_type}-{child_id}")
            if parent_id is None:
                raise MissingParentError(f"Can't find {parent_type} parent for {child_id}")
            return parent_id
        raise MissingParentError(f"Unsupported parent type {parent_type}")

    def add_gene(self, feature: SeqFeature) -> None:
        """Add the functional annotation of a given gene.

        Raises:
            DuplicateIdError: Do not allow genes with the same ID.

        """
        if feature.id in self.genes:
            raise DuplicateIdError(f"Gene ID {feature.id} already added")
        gene = self._generic_feature(feature, "gene")
        self.genes[gene["id"]] = gene

    def add_transcript(self, feature: SeqFeature, gene_id: str) -> None:
        """Add the functional annotation of a given transcript,
        and keep a record of its gene parent.

        Raises:
            DuplicateIdError: Do not allow transcripts with the same ID.

        """
        if feature.id in self.transcripts:
            raise DuplicateIdError(f"Transcript ID {feature.id} already added")
        transcript = self._generic_feature(feature, "transcript")
        self.transcripts[transcript["id"]] = transcript
        self.add_parent("gene", gene_id, transcript["id"])

    def add_translation(self, feature: SeqFeature, transcript_id: str) -> None:
        """Add the functional annotation of a given translation,
        and keep a record of its transcript parent.

        Raises:
            DuplicateIdError: Do not allow translations with the same ID.

        """
        if feature.id in self.translations:
            raise DuplicateIdError(f"Translation ID {feature.id} already added")
        translation = self._generic_feature(feature, "translation")
        self.translations[translation["id"]] = translation
        self.parent[f"transcript-{translation['id']}"] = transcript_id
        self.add_parent("transcript", transcript_id, translation["id"])

    def add_transposable_element(self, feature: SeqFeature) -> None:
        """Add the functional annotation of a transposable element.

        Raises:
            DuplicateIdError: Do not allow transposable elements with the same ID.

        """
        if feature.id in self.transposable_elements:
            raise DuplicateIdError(f"Transposable element ID {feature.id} already added")
        te = self._generic_feature(feature, "transposable_element")
        self.transposable_elements[te["id"]] = te

    def _generic_feature(self, feature: SeqFeature, feat_type: str) -> Dict[str, Any]:
        """Create a feature object following the specifications.

        Args:
            feature: The SeqFeature to add to the list.
            feat_type: Feature type of the feature to store (e.g. gene, transcript, translation).

        """

        feature_object: Annotation = {"object_type": feat_type, "id": feature.id}

        # Description?
        if "product" in feature.qualifiers:
            description = feature.qualifiers["product"][0]
            if self.product_is_informative(description):
                feature_object["description"] = description

        if "Name" in feature.qualifiers and "description" not in feature_object:
            feature_object["description"] = feature.qualifiers["Name"][0]

        # Don't keep useless description
        if ("description" in feature_object) and not self.product_is_informative(
            feature_object["description"], feature.id
        ):
            del feature_object["description"]

        # Synonyms?
        if "Name" in feature.qualifiers:
            feat_name = feature.qualifiers["Name"][0]
            if feat_name != feature.id:
                feature_object["synonyms"] = {"synonym": feat_name, "default": True}

        # is_pseudogene?
        if feature.type.startswith("pseudogen"):
            feature_object["is_pseudogene"] = True

        return feature_object

    def _transfer_descriptions(self) -> None:
        # Transfer translation CDS if useful
        for translation_id, translation in self.translations.items():
            description = translation.get("description")
            if description is not None:
                # Check transcript
                parent_tr_id = self.get_parent("transcript", translation_id)
                parent_tr = self.transcripts[parent_tr_id]
                tr_description = parent_tr.get("description")
                if tr_description is None:
                    parent_tr["description"] = description

        # Same, from transcripts to genes
        for transcript_id, transcript in self.transcripts.items():
            description = transcript.get("description")
            if description is not None:
                # Check gene
                parent_gene_id = self.get_parent("gene", transcript_id)
                parent_gene = self.genes[parent_gene_id]
                tr_description = parent_gene.get("description")
                if tr_description is None:
                    parent_gene["description"] = description

    @staticmethod
    def product_is_informative(product: str, feat_id: Optional[str] = None) -> bool:
        """Returns True if the product name contains informative words, False otherwise.

        It is considered uninformative when the description contains words such as "hypothetical" or
        or "putative". If a feature ID is provided, consider it uninformative as well (we do not want
        descriptions to be just the ID).

        Args:
            product: A product name.
            feat_id: Feature ID (optional).

        """
        non_informative_words = [
            "hypothetical",
            "putative",
            "uncharacterized",
            "unspecified",
            r"(of )?unknown function",
            "conserved",
            "predicted",
            "fragment",
            "product",
            "protein",
            "RNA",
            r"variant( \d+)?",
        ]
        non_informative_re = re.compile(r"|".join(non_informative_words), re.IGNORECASE)

        # Remove the feature ID if it's in the description
        if feat_id is not None:
            feat_id_re = re.compile(feat_id, re.IGNORECASE)
            product = re.sub(feat_id_re, "", product)

        # Remove punctuations
        punct_re = re.compile(r"[,;: _()-]+")
        product = re.sub(punct_re, " ", product)

        # Then remove non informative words
        product = re.sub(non_informative_re, " ", product)

        # Anything (informative) left?
        empty_re = re.compile(r"^[ ]*$")
        if empty_re.match(product):
            return False
        return True

    def _to_list(self):
        feat_list = []
        for feat in (
            list(self.genes.values())
            + list(self.transcripts.values())
            + list(self.translations.values())
            + list(self.transposable_elements.values())
        ):
            feat_list.append(feat)
        return feat_list

    def to_json(self, out_path: PathLike) -> None:
        """Print out the current annotation list in a json file.

        Args:
            out_path: JSON file path where to write the data.

        """
        self._transfer_descriptions()
        feats_list = self._to_list()
        print_json(Path(out_path), feats_list)
