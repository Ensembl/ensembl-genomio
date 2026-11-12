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
"""Methods for querying the Ensembl core database for analysis data."""

from sqlalchemy import select
from sqlalchemy.engine import URL

from ensembl.core.models import Analysis, AnalysisDescription

from ensembl.utils.database import DBConnection

__all__ = ["get_analysis_description"]


def get_analysis_description(
    logic_name: str,
    prod_url: URL,
) -> tuple[str, str]:
    """Retrieve analysis description metadata from the production database.

    Args:
        logic_name: Analysis logic name to look up.
        prod_url: URL for the production database.

    Returns:
        A tuple containing the analysis description and display label.

    Raises:
        ValueError: If the analysis description is not found in the production database.

    """
    database = DBConnection(prod_url, reflect=False)

    with database.session_scope() as session:
        result = session.execute(
            select(
                AnalysisDescription.description,
                AnalysisDescription.display_label,
            )
            .join(Analysis, Analysis.analysis_id == AnalysisDescription.analysis_id)
            .where(Analysis.logic_name == logic_name)
        ).one_or_none()

        if result is None:
            raise ValueError(f"Analysis not found in production database: {logic_name}")

        return result.description, result.display_label
