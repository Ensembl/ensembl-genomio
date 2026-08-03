"""Literature-based genome-assembly metadata extraction.

Given an NCBI assembly accession, find the assembly's publication, read its
text (and supplementary files), and extract species, ploidy, chromosome
number, cultivar/strain and sex via a rule + vector + optional-LLM ensemble.
"""
from .pipeline import run_pipeline, run_batch

__all__ = ["run_pipeline", "run_batch"]
