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
"""Parse publication text (full-text XML, abstract, supplementary) into weighted sections and chunks."""

import logging
import re
from xml.etree import ElementTree as ET

logger = logging.getLogger(__name__)

_TEX_BLOCK = re.compile(r"\\documentclass.*?\\end\{document\}", re.DOTALL)
_TEX_CMD = re.compile(r"\\[a-zA-Z]+\s*(\[[^\]]*\])?(\{[^}]*\})?")
_HTML_TAG = re.compile(r"<[^>]+>")
# Match ALL whitespace incl. Unicode (thin space U+2009, nbsp U+00A0, …). Genome
# papers often typeset chromosome formulas as "2n = 4x = 30"
# with thin spaces; normalising them to plain spaces lets the extractor's regex
# read the formula.
_WS = re.compile(r"\s+")


def clean_text(text: str) -> str:
    """Strip LaTeX/HTML markup and collapse all whitespace (including Unicode thin
    spaces / nbsp) to single spaces. Returns '' for empty input."""
    if not text:
        return ""
    text = _TEX_BLOCK.sub(" ", text)
    text = _TEX_CMD.sub(" ", text)
    text = text.replace("$$", " ").replace("$", " ")
    return _WS.sub(" ", text).strip()


def detect_text_type(text: str) -> str:
    """Classify raw text as full-text 'xml' (starts with a tag) or plain 'abstract'."""
    return "xml" if text.lstrip().startswith("<") else "abstract"


def _serialize_table(table_wrap: ET.Element) -> str:
    """Flatten a <table-wrap> element into a compact "header: value | ..." string
    (first 40 rows) so chromosome / karyotype tables become searchable text."""
    table = table_wrap.find(".//table")
    if table is None:
        return ""

    headers = []
    thead = table.find(".//thead")
    header_row = thead.find(".//tr") if thead is not None else None
    if header_row is None:
        header_row = table.find(".//tr")
    if header_row is not None:
        headers = [clean_text("".join(c.itertext())) for c in list(header_row) if c.tag in ("th", "td")]

    tbody = table.find(".//tbody")
    row_elems = tbody.findall(".//tr") if tbody is not None else table.findall(".//tr")

    rows = []
    for tr in row_elems:
        cells = [clean_text("".join(c.itertext())) for c in list(tr) if c.tag in ("th", "td")]
        if not any(cells):
            continue
        if headers and len(headers) == len(cells):
            row_text = " | ".join(f"{header}: {value}" for header, value in zip(headers, cells) if value)
        else:
            row_text = " | ".join(value for value in cells if value)
        if row_text:
            rows.append(row_text)

    return " ; ".join(rows[:40])


def parse_xml_tables(root: ET.Element) -> dict[str, str]:
    """Extract each <table-wrap> from a PMC article as a {title: caption + rows}
    section, keeping the longest body when a title repeats."""
    tables: dict[str, str] = {}
    for i, tw in enumerate(root.iter("table-wrap"), 1):
        label = clean_text(tw.findtext("label", ""))

        caption_texts = []
        for cap in tw.iter("caption"):
            for p in cap.iter("p"):
                t = clean_text("".join(p.itertext()))
                if t:
                    caption_texts.append(t)
        caption = " ".join(caption_texts)

        body = clean_text(f"{caption} {_serialize_table(tw)}")
        if not body:
            continue

        title = label or f"Table {i}"
        if title in tables and len(tables[title]) >= len(body):
            continue
        tables[title] = body

    return tables


def parse_xml_sections(xml_text: str) -> dict[str, str]:
    """Parse PMC full-text XML into a {section title: text} dict, including the
    abstract, body sections (longest kept on title clashes) and tables."""
    sections = {}
    try:
        root = ET.fromstring(xml_text)

        for abstract in root.iter("abstract"):
            texts = []
            for p in abstract.iter("p"):
                full_text = clean_text("".join(p.itertext()))
                if full_text:
                    texts.append(full_text)
            if texts:
                sections["Abstract"] = " ".join(texts)

        for sec in root.iter("sec"):
            title_elem = sec.find("title")
            if title_elem is None or not title_elem.text:
                continue
            section_title = title_elem.text.strip()
            texts = []
            for p in sec.iter("p"):
                full_text = clean_text("".join(p.itertext()))
                if full_text:
                    texts.append(full_text)
            if texts:
                existing = sections.get(section_title, "")
                merged = " ".join(texts)
                sections[section_title] = merged if len(merged) > len(existing) else existing

        for title, body in parse_xml_tables(root).items():
            sections[title] = body

    except ET.ParseError as e:
        logger.warning(f"  [parse] XML parse error: {e}")

    return sections


def parse_abstract(abstract_text: str) -> dict[str, str]:
    """Wrap a plain abstract string as a single {"Abstract": cleaned text} section."""
    # Abstract snippets from search results sometimes carry stray markup
    # (e.g. "<title>Abstract</title> <p>..."); strip tags before cleaning so
    # the section body is plain text.
    text = _HTML_TAG.sub(" ", abstract_text or "")
    return {"Abstract": clean_text(text)}


def parse_supplementary(supp_text: str) -> dict[str, str]:
    """Turn the concatenated supplementary blob produced by
    fetch.get_supplementary_text_by_pmcid into sections, one per file. Each block
    starts with a '[Supplementary file: <name>]' marker. Section titles are
    prefixed 'Supplementary:' so downstream section weighting can recognise them."""
    if not supp_text:
        return {}
    sections = {}
    blocks = re.split(r"\[Supplementary file:\s*([^\]]+)\]", supp_text)
    # re.split keeps capture groups: [pre, name1, body1, name2, body2, ...]
    for i in range(1, len(blocks) - 1, 2):
        name = blocks[i].strip()
        body = clean_text(blocks[i + 1])
        if body:
            sections[f"Supplementary: {name}"] = body
    if not sections:  # no markers found — keep as one section
        body = clean_text(supp_text)
        if body:
            sections["Supplementary"] = body
    return sections


def split_into_chunks(sections: dict[str, str], chunk_size: int = 300) -> list[dict]:
    """Split each section's text into sentence-aware chunks of ~chunk_size chars,
    returning a list of {section, text, chunk_index} dicts for embedding/search."""
    # chunk_size is in characters: ~300 keeps each chunk to a few sentences, which
    # suits the sentence-transformer in search.py (short, focused passages embed
    # and retrieve better than long ones) while staying well under its token limit.
    chunks = []
    for sec_title, text in sections.items():
        sentences = text.replace(". ", ".|").split("|")
        current_chunk = ""
        chunk_index = 0

        for sent in sentences:
            if len(current_chunk) + len(sent) < chunk_size:
                current_chunk += sent + " "
            else:
                if current_chunk.strip():
                    chunks.append(
                        {
                            "section": sec_title,
                            "text": current_chunk.strip(),
                            "chunk_index": chunk_index,
                        }
                    )
                    chunk_index += 1
                current_chunk = sent + " "

        if current_chunk.strip():
            chunks.append(
                {
                    "section": sec_title,
                    "text": current_chunk.strip(),
                    "chunk_index": chunk_index,
                }
            )

    return chunks


def parse_paper(text_data: dict) -> dict:
    """Parse one paper's text_data (full-text XML, abstract, or supplementary blob)
    into weighted sections and search chunks. Returns a dict with source, sections,
    chunks, n_sections and n_chunks."""
    source = text_data.get("source", "none")
    text = text_data.get("text")

    if not text:
        logger.info("  [parse] No text available — returning empty result")
        return {
            "source": source,
            "sections": {},
            "chunks": [],
            "n_sections": 0,
            "n_chunks": 0,
        }

    if detect_text_type(text) == "xml":
        sections = parse_xml_sections(text)
        logger.info(f"  [parse] XML parsed → {len(sections)} sections")
        # Fallback: an abstract snippet that merely *starts* with a tag
        # (e.g. "<title>Abstract</title> <p>...") is not a full PMC article, so
        # the XML parser yields nothing. Treat it as an abstract instead.
        if not sections:
            sections = parse_abstract(text)
            logger.info(
                f"  [parse] XML yielded 0 sections → fell back to abstract " f"({len(sections)} section)"
            )
    else:
        sections = parse_abstract(text)
        logger.info(f"  [parse] Plain abstract wrapped as single section")

    # Fold in supplementary-file text (chromosome/karyotype tables often live
    # here, not in the article body). Attached by fetch.get_best_available_text.
    supp_sections = parse_supplementary(text_data.get("supplementary", ""))
    if supp_sections:
        sections.update(supp_sections)
        logger.info(f"  [parse] + {len(supp_sections)} supplementary section(s)")

    chunks = split_into_chunks(sections)
    logger.info(f"  [parse] {len(chunks)} chunks ready for extract.py and search.py")

    return {
        "source": source,
        "sections": sections,
        "chunks": chunks,
        "n_sections": len(sections),
        "n_chunks": len(chunks),
    }
