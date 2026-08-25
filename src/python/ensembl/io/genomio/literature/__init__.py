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
"""Literature-based genome-assembly metadata extraction.

Given an NCBI assembly accession, find the assembly's publication, read its
text (and supplementary files), and extract species, ploidy, chromosome
number, cultivar/strain and sex via a rule + vector + optional-LLM ensemble.
"""

__all__ = ["run_batch", "run_pipeline"]

from ensembl.io.genomio.literature.pipeline import run_batch, run_pipeline
