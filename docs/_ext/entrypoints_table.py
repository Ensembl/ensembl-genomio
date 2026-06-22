# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Sphinx extension: auto-inject a CLI entry-points reference table into a Markdown page.

At ``builder-inited`` time this extension:

1. Reads ``[project.scripts]`` from the ``pyproject.toml``.
2. Builds a Markdown table mapping each script name to its Python target (``package.module:function``).
3. Writes the table between two sentinel comments inside the target Markdown file, replacing whatever
   was there from a previous build so the file stays under version control and is always up-to-date.

Usage:
    Add the extension directory to ``sys.path`` in ``conf.py`` and include the extension::

        import sys, pathlib

        sys.path.insert(0, str(pathlib.Path(__file__).parent / "_ext"))
        extensions = [..., "entrypoints_table"]

    Optionally override the defaults in ``conf.py``::

        entrypoints_target_file = "user_guide/usage.md"  # path relative to docs source dir
        entrypoints_toml_file = "../pyproject.toml"  # path relative to docs source dir

"""

import logging
import re
import sys
from pathlib import Path
from typing import Any

from sphinx.application import Sphinx
from sphinx.util import logging as sphinx_logging

if sys.version_info >= (3, 11):
    import tomllib
try:
    import tomli as tomllib
except ImportError as exc:
    raise ImportError("Python < 3.11 requires the 'tomli' package: pip install tomli") from exc


logger: logging.Logger = sphinx_logging.getLogger(__name__)

# Sentinels written into the Markdown file. Everything between them is replaced.
_SENTINEL_START = "<!-- entrypoints-table:start -->"
_SENTINEL_END = "<!-- entrypoints-table:end -->"

# Regex that matches the region between the two sentinels (inclusive).
_REGION_RE = re.compile(rf"{re.escape(_SENTINEL_START)}.*?{re.escape(_SENTINEL_END)}", re.DOTALL)


def _read_entry_points(toml_path: Path) -> dict[str, str]:
    """Parse ``project.scripts`` from a ``pyproject.toml`` file.

    Args:
        toml_path: Absolute path to the ``pyproject.toml`` file.

    Returns:
        Mapping of ``script-name`` to ``package.module:function``. Returns an empty dict when
        the section is absent.

    """
    with toml_path.open("rb") as fh:
        data = tomllib.load(fh)
    return data.get("project", {}).get("scripts", {})


def _build_markdown_table(entry_points: dict[str, str]) -> str:
    """Render ``entry_points`` as a Markdown table string.

    The table has two columns: **Command** and **Python target**.

    Args:
        entry_points: Mapping returned by :func:`_read_entry_points`.

    Returns:
        A complete Markdown table (no trailing newline), or a short italicised notice when
        ``entry_points`` is empty.

    """
    if not entry_points:
        return "_No entry points are defined in `pyproject.toml`._"
    # Column widths for pretty-printing (purely cosmetic for raw Markdown readers).
    max_cmd = max(len(cmd) for cmd in entry_points)
    max_target = max(len(target) for target in entry_points.values())
    col_cmd = max(len("Command"), max_cmd)
    col_target = max(len("Python target"), max_target)
    header = f"| {'Command':<{col_cmd}} | {'Python target':<{col_target}} |"
    separator = f"| {'-' * col_cmd} | {'-' * col_target} |"
    rows = [
        f"| `{cmd}`{' ' * (col_cmd - len(cmd) - 2)} | `{target}`{' ' * (col_target - len(target) - 2)} |"
        for cmd, target in sorted(entry_points.items())
    ]
    return "\n".join([header, separator, *rows])


def _inject_table(usage_path: Path, table_md: str) -> None:
    """Replace the sentinel region in *usage_path* with *table_md*.

    If the sentinels are not found the table (including the sentinels) is appended to the end of the file
    with a preceding blank line, so the extension is safe to add to a page that has not yet been prepared.

    Args:
        usage_path: Absolute path to the target Markdown file.
        table_md: Rendered Markdown table produced by :func:`_build_markdown_table`.

    """
    original = usage_path.read_text(encoding="utf-8")
    replacement_block = f"{_SENTINEL_START}\n{table_md}\n{_SENTINEL_END}"
    if _REGION_RE.search(original):
        updated = _REGION_RE.sub(replacement_block, original)
    else:
        logger.warning(
            "entrypoints_table: sentinels not found in '%s'. Appending the table at the end of the file.",
            usage_path,
        )
        updated = original.rstrip("\n") + f"\n\n{replacement_block}\n"
    usage_path.write_text(updated, encoding="utf-8")


def _on_builder_inited(app: Sphinx) -> None:
    """Sphinx event handler for ``builder-inited``.

    Reads configuration values set in ``conf.py``, resolves file paths relative to the Sphinx source
    directory, and orchestrates the table generation and injection.

    Args:
        app: The Sphinx application object provided by the event system.

    """
    src_dir = Path(app.srcdir)
    target_rel: str = getattr(app.config, "entrypoints_target_file", "usage.md")
    toml_rel: str = getattr(app.config, "entrypoints_toml_file", "../pyproject.toml")
    usage_path = src_dir / target_rel
    toml_path = (src_dir / toml_rel).resolve()
    logger.info("entrypoints_table: reading entry points from '%s'", toml_path)
    try:
        entry_points = _read_entry_points(toml_path)
    except FileNotFoundError as exc:
        logger.warning("entrypoints_table: %s — skipping table generation.", exc)
        return
    table_md = _build_markdown_table(entry_points)
    if not usage_path.is_file():
        logger.warning(
            "entrypoints_table: target file '%s' does not exist — skipping injection.",
            usage_path,
        )
        return
    _inject_table(usage_path, table_md)
    logger.info("entrypoints_table: table injected into '%s'.", usage_path)


def setup(app: Sphinx) -> dict[str, Any]:
    """Register the extension with Sphinx.

    Adds two optional ``conf.py`` configuration values:

    - ``entrypoints_target_file``: Path to the Markdown file to inject into, relative to the Sphinx
      source directory.
    - ``entrypoints_toml_file``: Path to ``pyproject.toml``, relative to the Sphinx source directory.

    Args:
        app: The Sphinx application object.

    Returns:
        Sphinx extension metadata.

    """
    app.add_config_value("entrypoints_target_file", "user_guide/usage.md", "env")
    app.add_config_value("entrypoints_toml_file", "../pyproject.toml", "env")
    app.connect("builder-inited", _on_builder_inited)
    return {
        "version": "1.0.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
