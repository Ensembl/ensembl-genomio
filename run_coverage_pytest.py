#!/usr/bin/env python3
"""
Run pytest under coverage, but pre-import Biopython before coverage starts.

This avoids a Biopython/NumPy optional-dependency check failing only when pytest is
bootstrapped under coverage tracing in some environments.
"""

from __future__ import annotations

import sys

import coverage


def main(argv: list[str]) -> int:
    # 1) Warm Biopython import *before* coverage starts
    #    (conftest imports Bio.SeqIO at import-time)
    from Bio import SeqIO  # noqa: F401

    # 2) Start coverage (adjust `source=` to your package root)
    cov = coverage.Coverage(
        branch=True,
        source=["ensembl"],  # or ["ensembl.io.genomio"] if you prefer narrower
    )
    cov.start()

    # 3) Run pytest
    import pytest

    exit_code = pytest.main(argv)

    # 4) Stop + save coverage
    cov.stop()
    cov.save()

    return exit_code


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
