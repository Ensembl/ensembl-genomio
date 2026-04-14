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
"""Unit testing of `ensembl.io.genomio.assembly.status` module."""
# pylint: disable=too-many-positional-arguments

from __future__ import annotations

from contextlib import nullcontext as does_not_raise
import json
import os
from pathlib import Path
from typing import Any, Callable, ContextManager
from unittest.mock import Mock, patch

from deepdiff import DeepDiff
import pytest
from pytest import FixtureRequest, param, raises
from sqlalchemy import Column, Index, String, text
from sqlalchemy.dialects.mysql import INTEGER
from sqlalchemy.engine import make_url
from sqlalchemy.orm import declarative_base

from ensembl.io.genomio.assembly.status import (
    DATASETS_SINGULARITY,
    UnsupportedFormatError,
    ReportStructure,
    extract_assembly_metadata,
    fetch_datasets_reports,
    fetch_accessions_from_core_dbs,
    generate_report_tsv,
    get_assembly_accessions,
    singularity_image_setter,
)
from ensembl.io.genomio.utils.json_utils import get_json
from ensembl.utils import StrPath
from ensembl.utils.database import UnitTestDB


MINIMUM_METADATA = {
    "my_core": ReportStructure(
        species_name="Octopus bimaculoides",
        taxon_id=37653,
        assembly_name="ASM119413v2",
        assembly_type="haploid",
        accession="GCF_001194135.2",
        assembly_status="current",
    )
}

STRAIN_METADATA = {
    "my_core": ReportStructure(
        species_name="Octopus bimaculoides",
        taxon_id=37653,
        strain="UCB-OBI-ISO-001",
        assembly_name="ASM119413v2",
        assembly_type="haploid",
        accession="GCF_001194135.2",
        assembly_status="current",
    )
}

COMPLETE_METADATA = {
    "my_core": ReportStructure(
        species_name="Octopus bimaculoides",
        taxon_id=37653,
        strain="UCB-OBI-ISO-001",
        assembly_name="ASM119413v2",
        assembly_type="haploid",
        accession="GCF_001194135.2",
        paired_assembly="GCA_001194135.2",
        last_updated="2015-06-29T09:51:41.073",
        assembly_status="current",
        assembly_notes="RefSeq",
    )
}


Base: Any = declarative_base()  # only possible type hint compatible with SQLAlchemy 1.4 and 2.0+


class Meta(Base):
    """Meta class mirroring the Ensembl core database meta table without any foreign keys"""

    __tablename__ = "meta"
    __table_args__ = (
        Index("species_value_idx", "species_id", "meta_value"),
        Index("species_key_value_idx", "species_id", "meta_key", "meta_value", unique=True),
    )

    meta_id: Column = Column(INTEGER(11), primary_key=True)
    species_id: Column = Column(INTEGER(10), nullable=False, index=True, server_default=text("'1'"))
    meta_key: Column = Column(String(40), nullable=False)
    meta_value: Column = Column(String(255), nullable=False)


@pytest.mark.dependency()
def test_report_structure() -> None:
    """Tests the `ReportStructure` class."""
    assert ReportStructure()


@pytest.mark.dependency(depends=["test_report_structure"])
def test_report_structure_to_dict() -> None:
    """Tests the `ReportStructure.to_dict()` method."""
    assert not DeepDiff(
        COMPLETE_METADATA["my_core"].to_dict(),
        {
            "Species Name": "Octopus bimaculoides",
            "Taxon ID": "37653",
            "Isolate/Strain": "UCB-OBI-ISO-001",
            "Asm name": "ASM119413v2",
            "Assembly type": "haploid",
            "Asm accession": "GCF_001194135.2",
            "Paired assembly": "GCA_001194135.2",
            "Asm last updated": "2015-06-29T09:51:41.073",
            "Asm status": "current",
            "Asm notes": "RefSeq",
        },
    )


@pytest.mark.dependency(depends=["test_report_structure"])
def test_report_structure_header() -> None:
    """Tests the `ReportStructure.header()` method."""
    expected_header = [
        "Species Name",
        "Taxon ID",
        "Isolate/Strain",
        "Asm name",
        "Assembly type",
        "Asm accession",
        "Paired assembly",
        "Asm last updated",
        "Asm status",
        "Asm notes",
    ]
    assert COMPLETE_METADATA["my_core"].header() == expected_header


@pytest.mark.dependency(depends=["test_report_structure"])
def test_report_structure_values() -> None:
    """Tests the `ReportStructure.values()` method."""
    expected_values = [
        "Octopus bimaculoides",
        "37653",
        "UCB-OBI-ISO-001",
        "ASM119413v2",
        "haploid",
        "GCF_001194135.2",
        "GCA_001194135.2",
        "2015-06-29T09:51:41.073",
        "current",
        "RefSeq",
    ]
    assert COMPLETE_METADATA["my_core"].values() == expected_values


@pytest.mark.parametrize(
    "sif_cache_dir, datasets_version, nextflow_cachedir, singularity_cachedir",
    [
        param("sif_cache", None, None, None, id="Personal SIF cache, default datasets"),
        param(None, None, "nxf_cache", None, id="Nextflow SIF cache, default datasets"),
        param(None, None, None, "singularity_cache", id="Singularity SIF cache, default datasets"),
        param(None, "my_datasets", None, None, id="No SIF cache, user datasets"),
    ],
)
@patch("ensembl.io.genomio.assembly.status.Client")
def test_singularity_image_setter(
    mock_client: Mock,
    tmp_path: Path,
    sif_cache_dir: Path | None,
    datasets_version: str | None,
    nextflow_cachedir: str | None,
    singularity_cachedir: str | None,
) -> None:
    """Tests the `singularity_image_setter()` function.

    Fixtures: tmp_path

    Args:
        sif_cache_dir: Path to locate existing, or download new SIF container image.
        datasets_version: URL of singularity container (custom `datasets` version if desired).
        nextflow_cachedir: Value to assign to environment variable NXF_SINGULARITY_CACHEDIR.
        singularity_cachedir: Value to assign to environment variable SINGULARITY_CACHEDIR.
    """
    mock_client.pull.return_value = True
    # Define SIF cache path and expected path used to pull the container
    sif_cache_path = None
    if sif_cache_dir:
        sif_cache_path = tmp_path / sif_cache_dir
        sif_cache_path.mkdir()
        expected_cache_path = sif_cache_path
    elif nextflow_cachedir:
        expected_cache_path = Path(nextflow_cachedir)
    elif singularity_cachedir:
        expected_cache_path = Path(singularity_cachedir)
    else:
        expected_cache_path = Path()
    # Get expected container URL used to pull the container
    if datasets_version:
        expected_container_url = datasets_version
    else:
        expected_container_url = DATASETS_SINGULARITY["datasets_version_url"]
    # Patch the environment variables
    new_env = {}
    if nextflow_cachedir:
        new_env["NXF_SINGULARITY_CACHEDIR"] = nextflow_cachedir
    if singularity_cachedir:
        new_env["SINGULARITY_CACHEDIR"] = singularity_cachedir
    with patch.dict(os.environ, new_env, clear=True):
        assert singularity_image_setter(sif_cache_path, datasets_version)
    # Check that the spython pull method was called with the right arguments
    mock_client.pull.assert_called_with(
        expected_container_url, stream=False, pull_folder=expected_cache_path, quiet=True
    )


@pytest.mark.parametrize(
    "file_name, expected_output, expectation",
    [
        param("assemblies.txt", ["GCA_900524668.2", "GCF_123456789.1"], does_not_raise(), id="Valid file"),
        param("meta_report.tsv", [], raises(UnsupportedFormatError), id="Wrong accession format"),
    ],
)
def test_get_assembly_accessions(
    data_dir: Path, file_name: str, expected_output: list[str], expectation: ContextManager
) -> None:
    """Tests the `get_assembly_accessions()` function.

    Fixtures:
        data_dir

    Args:
        file_name: File with one line per INSDC assembly accession.
        expected_output: Expected assembly accessions returned.
        expectation: Context manager of expected raise exception.
    """
    file_path = data_dir / file_name
    with expectation:
        result = get_assembly_accessions(file_path)
        assert result == expected_output


@pytest.mark.parametrize(
    "test_dbs",
    [
        [
            {"src": "one_core_db", "metadata": Base.metadata},
            {"src": "another_core_db", "metadata": Base.metadata},
            {"src": "yet_another_core_db", "metadata": Base.metadata},
        ],
    ],
    indirect=True,
)
def test_fetch_accessions_from_core_dbs(
    request: FixtureRequest, tmp_path: Path, test_dbs: dict[str, UnitTestDB]
) -> None:
    """Tests the `fetch_accessions_from_core_dbs()` function.

    Fixtures: request, tmp_path, test_dbs
    """
    # Create a file with each test database name
    tmp_file = tmp_path / "core_dbs_list.txt"
    with tmp_file.open("w") as fin:
        for db in test_dbs.values():
            fin.write(f"{db.dbc.db_name}\n")
    test_server = make_url(request.config.getoption("server"))
    accessions = fetch_accessions_from_core_dbs(tmp_file, test_server)
    assert not DeepDiff(accessions, {test_dbs["one_core_db"].dbc.db_name: "GCF_001194135.2"})


@patch("ensembl.io.genomio.assembly.status.Client")
def test_fetch_datasets_reports(
    mock_client: Mock, tmp_path: Path, data_dir: Path, assert_files: Callable[[StrPath, StrPath], None]
) -> None:
    """Tests the `fetch_datasets_reports()` function.

    Fixtures:
        tmp_path, data_dir, assert_files
    """

    def execute_return(
        command: list[str], **kwargs: Any  # pylint: disable=unused-argument
    ) -> dict[str, str]:
        report_path = data_dir / f"{command[-1]}.asm_report.json"
        if report_path.exists():
            report = get_json(report_path)
            return {"message": f'{{"reports": [{json.dumps(report)}], "total_count": 1}}'}
        return {"message": '{"total_count": 0}'}

    mock_client.execute.side_effect = execute_return
    accessions = {"my_core": "GCF_001194135.2", "your_core": "GCA_001122334.1"}
    report = fetch_datasets_reports(mock_client, accessions, tmp_path, 1)
    assert not (tmp_path / "GCA_001122334.1.asm_report.json").exists()
    report_file = tmp_path / "GCF_001194135.2.asm_report.json"
    expected_report_file = data_dir / "GCF_001194135.2.asm_report.json"
    assert_files(report_file, expected_report_file)
    expected_report = get_json(expected_report_file)
    assert not DeepDiff(report, {"my_core": expected_report})


@patch("ensembl.io.genomio.assembly.status.Client")
def test_fetch_datasets_reports_value_error(mock_client: Mock) -> None:
    """Tests the `fetch_datasets_reports()` function when `ValueError` is raised."""
    mock_client.execute.return_value = {"message": [["unexpected nested list"]]}
    accessions = {"my_core": "GCF_001194135.2"}
    with raises(ValueError):
        fetch_datasets_reports(mock_client, accessions, Path(), 1)


@patch("ensembl.io.genomio.assembly.status.Client")
def test_fetch_datasets_reports_runtime_error(mock_client: Mock) -> None:
    """Tests the `fetch_datasets_reports()` function when `RuntimeError` is raised."""
    mock_client.execute.return_value = {"message": "FATAL error message"}
    accessions = {"my_core": "GCF_001194135.2"}
    with raises(RuntimeError, match=r"Singularity image execution failed! -> '.*'"):
        fetch_datasets_reports(mock_client, accessions, Path(), 1)


@pytest.mark.dependency(depends=["test_report_structure"])
@pytest.mark.parametrize(
    "file_name, expected_metadata",
    [
        param("simple.asm_report.json", MINIMUM_METADATA, id="Minimum report"),
        param("strain.asm_report.json", STRAIN_METADATA, id="Strain report"),
        param("cultivar.asm_report.json", MINIMUM_METADATA, id="Unexpected infraspecific name (cultivar)"),
        param("GCF_001194135.2.asm_report.json", COMPLETE_METADATA, id="Detailed report"),
    ],
)
def test_extract_assembly_metadata(
    data_dir: Path, file_name: str, expected_metadata: dict[str, ReportStructure]
) -> None:
    """Tests the `extract_assembly_metadata()` function.

    Fixtures: data_dir

    Args:
        file_name: Test data file to extract the assembly metadata from.
        expected_metadata: Expected key value pairs of source name <> assembly report.
    """
    report_path = data_dir / file_name
    report = {"my_core": get_json(report_path)}
    metadata = extract_assembly_metadata(report)
    assert not DeepDiff(metadata, expected_metadata)


@pytest.mark.dependency(depends=["test_report_structure"])
def test_generate_report_tsv(
    tmp_path: Path, data_dir: Path, assert_files: Callable[[StrPath, StrPath], None]
) -> None:
    """Tests the `generate_report_tsv()` function.

    Fixtures: tmp_path, data_dir, assert_files
    """
    generate_report_tsv(COMPLETE_METADATA, "core_db", tmp_path, "test_output")
    assert_files(tmp_path / "test_output.tsv", data_dir / "meta_report.tsv")
