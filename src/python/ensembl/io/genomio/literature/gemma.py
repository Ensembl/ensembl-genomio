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
"""Optional local Gemma-3 LLM layer for grounded ploidy, sex and strain/cultivar extraction."""

import json
import logging
import os
import re
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

# ---- controlled vocabulary (shared with extract.py) ----
try:
    from ensembl.io.genomio.literature.extract import _LEVEL_LABEL
except Exception:
    _LEVEL_LABEL = {
        1: "monoploid",
        2: "diploid",
        3: "triploid",
        4: "tetraploid",
        5: "pentaploid",
        6: "hexaploid",
        8: "octoploid",
    }

VALID_MECHANISMS = {"autopolyploid", "allopolyploid", "amphidiploid", "segmental_allopolyploid"}

# Sex vocabulary used to normalise the model's `sex` field.
_SEX_VALID = {"male", "female", "hermaphrodite", "monoecious", "dioecious", "not_applicable", "unknown"}
_SEX_MAP = {"hermaphroditic": "hermaphrodite"}

# Ploidy word (incl. allo-/auto- prefixed forms) -> integer level.
_WORD_LEVEL = {
    "haploid": 1,
    "monoploid": 1,
    "diploid": 2,
    "triploid": 3,
    "tetraploid": 4,
    "pentaploid": 5,
    "hexaploid": 6,
    "allohexaploid": 6,
    "autohexaploid": 6,
    "heptaploid": 7,
    "octoploid": 8,
    "octaploid": 8,
    "allooctoploid": 8,
    "decaploid": 10,
    "dodecaploid": 12,
}

# Mechanism keyword -> controlled term. Ordered most-specific first so a scan
# stops at the strongest match (e.g. "segmental allopolyploid" before "allopolyploid").
_MECH_KEYWORDS = [
    ("segmental allopolyploid", "segmental_allopolyploid"),
    ("segmental_allopolyploid", "segmental_allopolyploid"),
    ("amphidiploid", "amphidiploid"),
    ("amphitetraploid", "amphidiploid"),
    ("allotetraploid", "allopolyploid"),
    ("allohexaploid", "allopolyploid"),
    ("allooctoploid", "allopolyploid"),
    ("allopolyploid", "allopolyploid"),
    ("allopolyploidization", "allopolyploid"),
    ("autotetraploid", "autopolyploid"),
    ("autohexaploid", "autopolyploid"),
    ("autopolyploid", "autopolyploid"),
]

# Model confidence string -> base float (inferred results get an extra multiplier).
_CONF_MAP = {"high": 0.78, "medium": 0.62, "low": 0.45}

# ---- config ----
_MODELS_DIR = Path(__file__).resolve().parent / "models"
_DEFAULT_MODEL = str(_MODELS_DIR / "google_gemma-3-4b-it-Q4_K_M.gguf")
BASE_URL = os.environ.get("GEMMA_BASE_URL", "http://localhost:8000/v1")
MODEL = os.environ.get("GEMMA_MODEL", _DEFAULT_MODEL)
TIMEOUT = int(os.environ.get("GEMMA_TIMEOUT", "300"))
# Cap on passage text sent to the model. Set high enough to hold a WHOLE trusted
# paper (not just the retrieved vector chunks) — a typical full text is well under
# this, so in practice nothing is truncated; ~100k chars fits inside the default
# context below (a paper is read one accession at a time, so latency is acceptable).
# Lower it for speed on the GGUF backend if needed.
MAX_CHARS = int(os.environ.get("GEMMA_MAX_CHARS", "100000"))
# GGUF context window. Gemma-3-4B supports up to 128k; the default holds a full
# paper plus the system prompt and few-shot examples.
N_CTX = int(os.environ.get("GEMMA_N_CTX", "32768"))
N_THREADS = int(os.environ.get("GEMMA_N_THREADS", "32"))  # one NUMA node

# ---- GGUF backend (llama-cpp-python, quantized, fastest CPU option) ----
_gguf_model = None

# ---- HuggingFace direct backend (float32 fallback, ~16GB RAM) ----
_HF_MODEL_DIR = os.environ.get("GEMMA_HF_MODEL_DIR", str(_MODELS_DIR / "gemma-3-4b-it"))
_direct_model = None
_direct_tokenizer = None


def _gguf_available() -> bool:
    return MODEL.endswith(".gguf") and Path(MODEL).is_file()


def _weights_available() -> bool:
    return (Path(_HF_MODEL_DIR) / "config.json").is_file()


def _server_is_up() -> bool:
    try:
        requests.get(f"{BASE_URL}/models", timeout=3)
        return True
    except Exception:
        return False


def is_enabled() -> bool:
    """Return whether the optional Gemma layer should run. `GEMMA_ENABLED=1`/`0`
    forces it on/off; otherwise it auto-enables when a GGUF file, local HF weights,
    or a running inference server is detected."""
    env_gemma = os.environ.get("GEMMA_ENABLED")

    if env_gemma == "1":
        return True
    if env_gemma == "0":
        return False
    # Auto-detect: enabled if GGUF file, HF weights, or HTTP server available
    return _gguf_available() or _weights_available() or _server_is_up()


def _load_gguf_backend() -> None:
    """Load the GGUF (llama-cpp) model into the module singleton on first use."""
    global _gguf_model

    if _gguf_model is not None:
        return
    from llama_cpp import Llama

    logger.info(f"  [gemma-gguf] loading {Path(MODEL).name} (first call, may take ~10s) ...")
    _gguf_model = Llama(
        model_path=MODEL,
        n_ctx=N_CTX,
        n_threads=N_THREADS,
        chat_format="gemma",
        verbose=False,
    )
    logger.info("  [gemma-gguf] model ready")


def _call_gguf(species: str, passages: list, accession: str | None = None) -> dict | None:
    try:
        _load_gguf_backend()
    except Exception as e:
        logger.warning(f"  [gemma-gguf] load failed ({e}); will try next backend")
        return None
    messages = _build_messages(species, passages, accession=accession)

    if _gguf_model is None:
        logger.warning("  [gemma-gguf] model failed to initialize")
        return None

    try:
        resp = _gguf_model.create_chat_completion(
            messages=messages,
            temperature=0,
            max_tokens=600,
            response_format={"type": "json_object"},
        )
        content = resp["choices"][0]["message"]["content"]  # type: ignore[index]
    except Exception as e:
        logger.warning(f"  [gemma-gguf] inference failed ({e})")
        return None
    return _parse_json(content)  # type: ignore[arg-type]


def _load_direct_backend() -> None:
    global _direct_model, _direct_tokenizer
    if _direct_model is not None:
        return

    model_name = Path(_HF_MODEL_DIR).name
    logger.info(f"  [gemma-direct] loading {model_name} (first call, may take a minute) ...")

    import torch
    from transformers import AutoTokenizer, AutoModelForCausalLM

    _direct_model = AutoModelForCausalLM.from_pretrained(_HF_MODEL_DIR, torch_dtype=torch.float32).eval()
    _direct_tokenizer = AutoTokenizer.from_pretrained(_HF_MODEL_DIR)
    logger.info("  [gemma-direct] model ready")


def _call_direct(species: str, passages: list, accession: str | None = None) -> dict | None:
    try:
        _load_direct_backend()
    except Exception as e:
        logger.warning(f"  [gemma-direct] load failed ({e}); will try HTTP server")
        return None

    import torch

    messages = _build_messages(species, passages, accession=accession)

    if _direct_model is None or _direct_tokenizer is None:
        logger.warning("  [gemma-direct] model failed to initialize")
        return None

    encoded = _direct_tokenizer.apply_chat_template(
        messages,
        add_generation_prompt=True,
        return_tensors="pt",
        tokenize=True,
    )
    if hasattr(encoded, "input_ids"):
        input_ids = encoded.input_ids
        attention_mask = encoded.attention_mask  # type: ignore[union-attr]
    else:
        input_ids = encoded
        attention_mask = torch.ones_like(input_ids)  # type: ignore[arg-type]

    eos_ids = [_direct_tokenizer.eos_token_id]
    eot = _direct_tokenizer.convert_tokens_to_ids("<end_of_turn>")
    if isinstance(eot, int):
        eos_ids.append(eot)

    try:
        with torch.no_grad():
            output = _direct_model.generate(
                input_ids,
                attention_mask=attention_mask,
                max_new_tokens=600,
                do_sample=False,
                pad_token_id=_direct_tokenizer.eos_token_id,
                eos_token_id=eos_ids,
            )
    except Exception as e:
        logger.warning(f"  [gemma-direct] generate failed ({e})")
        return None

    new_ids = output[0][input_ids.shape[1] :]
    text = _direct_tokenizer.decode(new_ids, skip_special_tokens=True)
    return _parse_json(text)  # type: ignore[arg-type]


# JSON schema used for vLLM guided decoding (forces well-formed output).
_PLOIDY_SCHEMA = {
    "type": "object",
    "properties": {
        "reasoning": {"type": "string"},
        "target_organism": {"type": ["string", "null"]},
        "ploidy_level": {"type": ["integer", "null"]},
        "ploidy_category": {"type": ["string", "null"]},
        "mechanism": {"type": ["string", "null"]},
        "chromosome_formula": {"type": ["string", "null"]},
        "evidence_quotes": {"type": "array", "items": {"type": ["string", "null"]}},
        "evidence_type": {"type": "string"},
        "stated_explicitly": {"type": "boolean"},
        "confidence": {"type": "string"},
        "ambiguity_note": {"type": ["string", "null"]},
        "sex": {"type": ["string", "null"]},
        "strain_cultivar": {"type": ["string", "null"]},
    },
    "required": [
        "reasoning",
        "ploidy_level",
        "ploidy_category",
        "mechanism",
        "evidence_quotes",
        "evidence_type",
        "stated_explicitly",
        "confidence",
        "sex",
        "strain_cultivar",
    ],
}

_SYSTEM = """\
You are a senior plant genomics scientist and database curator with deep \
expertise in polyploidy, chromosome biology, and plant genome evolution. You are \
building a high-quality genome-assembly database by reading passages from a \
scientific paper and determining the ploidy of a specified TARGET organism.

WHAT COUNTS AS PLOIDY EVIDENCE (you reason like an expert)
- Chromosome formulas encode it directly: the coefficient of x is the ploidy. \
"2n = 2x = 14" -> diploid (2); "2n = 4x = 28" -> tetraploid (4); \
"2n = 6x = 42" -> hexaploid (6).
- Polyploidy terminology carries a level even without a bare number: \
"allohexaploid", "amphidiploid", "autotetraploid", "polyploid with three subgenomes".
- Subgenome counts reveal it: bread wheat with A, B and D subgenomes (from three \
diploid progenitors) is hexaploid; cotton with A and D subgenomes is allotetraploid.
- Duplication / hybridization history: "arose from hybridization of two diploid \
species" implies tetraploid; a "whole-genome duplication" doubles the base level.
- Genome-size comparisons are indicative but weaker: "genome ~3x that of its \
diploid relative" suggests raised ploidy — treat as LOW confidence unless \
corroborated.

PLOIDY LEVEL = NUMBER OF CHROMOSOME SETS (the x in the nx formula)
haploid 1, diploid 2, triploid 3, tetraploid 4, pentaploid 5, hexaploid 6, \
octoploid 8. Report the integer in `ploidy_level` and the word in `ploidy_category`. \
IMPORTANT: for "2n = 4x = 30", ploidy_level = 4 (the x value), NOT 30 (the 2n \
chromosome count). ploidy_level is always a small integer (1–16).

WHICH ORGANISM (critical)
Determine the ploidy of the TARGET organism only (given in the user message). \
Do NOT report the ploidy of its progenitors, ancestors, or sister species. \
A single paper often names organisms at different ploidies — pick the target. \
If no target is specified, treat the organism whose genome/assembly the paper \
presents as the target.

GROUNDING (critical)
Base your answer ONLY on evidence in the passages. Do NOT answer from memorized \
species->ploidy facts. Passage-grounded inference (reading a formula, counting \
subgenomes, using stated progenitor counts) is expected and encouraged; memorized \
facts are not. If the passages contain no ploidy-relevant evidence, return \
ploidy_level = null, even if you believe you know the species' ploidy.

SEX
Report the sex of the sequenced individual if the paper states it. \
Use one of: male, female, hermaphrodite, monoecious, dioecious, not_applicable, unknown. \
Use not_applicable for species with no distinct sexes; unknown if not mentioned.

STRAIN / CULTIVAR
Report the specific cultivar, strain, or accession name used for the assembly \
(e.g. "Chinese Spring", "TM-1", "Col-0"); null if not specified.

OUTPUT — a single JSON object and nothing else
Fill every field. Use JSON null (not the string "null") for unknowns. Put all \
reasoning inside the `reasoning` field so there is no text outside the JSON.

{
  "reasoning": "<brief, concrete: what ploidy-relevant evidence is present, which organism it refers to, whether that is the target, and how you derived the level>",
  "target_organism": "<organism name as it appears in the passages, or null>",
  "ploidy_level": <integer x in "2n=nx=..." (e.g. 4 for "2n=4x=30"), NOT the 2n count; or null>,
  "ploidy_category": "<one of: haploid, diploid, triploid, tetraploid, pentaploid, hexaploid, heptaploid, octoploid, polyploid, aneuploid; or null>",
  "mechanism": "<one of: autopolyploid, allopolyploid, amphidiploid, segmental_allopolyploid; or null>",
  "chromosome_formula": "<verbatim formula if present, e.g. '2n = 6x = 42'; else null>",
  "evidence_quotes": ["<most decisive sentence verbatim>", "<2nd supporting sentence verbatim>", "<3rd supporting sentence verbatim>"],
  "evidence_type": "<one of: explicit_statement, chromosome_formula, subgenome_count, progenitor_inference, genome_size_comparison, terminology, none>",
  "stated_explicitly": <true only if the text directly names the ploidy (e.g. "hexaploid", "6x"); false if inferred; false when ploidy_level is null>,
  "confidence": "<high, medium, or low>",
  "ambiguity_note": "<note any conflicting ploidy signals or uncertainty about target vs. relative; else null>",
  "sex": "<male, female, hermaphrodite, monoecious, dioecious, not_applicable, or unknown>",
  "strain_cultivar": "<cultivar/strain/accession name, or null>"
}

FLAG GUIDANCE
- evidence_quotes: list up to 3 verbatim sentences, most informative first. \
Use an empty list [] when ploidy_level is null.
- stated_explicitly: true for direct naming ("hexaploid", "6x"); false for expert \
inference (formula math, subgenome/progenitor counting, size comparison).
- confidence: high = explicit statement or unambiguous chromosome formula for the \
target; medium = solid inference (subgenome or progenitor counts, clear terminology); \
low = indirect signals (genome-size ratios), partial evidence, or ploidy_level = null.
- mechanism: use null for diploids or when hybridization mode is unclear."""


_FEWSHOT = [
    {
        "target": "Triticum aestivum",
        "passage": (
            "Bread wheat (Triticum aestivum) is an allohexaploid (2n = 6x = 42) "
            "that arose from hybridization among three diploid progenitors "
            "carrying the A, B and D genomes. Here we present a chromosome-scale "
            "assembly of T. aestivum cv. Chinese Spring."
        ),
        "output": {
            "reasoning": (
                "The passage states T. aestivum is an allohexaploid and gives "
                "2n = 6x = 42 (x-coefficient 6 -> hexaploid). The A/B/D progenitors "
                "are diploid relatives, not the target. The assembly is of T. aestivum cv. Chinese Spring."
            ),
            "target_organism": "Triticum aestivum",
            "ploidy_level": 6,
            "ploidy_category": "hexaploid",
            "mechanism": "allopolyploid",
            "chromosome_formula": "2n = 6x = 42",
            "evidence_quotes": [
                "Bread wheat (Triticum aestivum) is an allohexaploid (2n = 6x = 42) that arose from hybridization among three diploid progenitors carrying the A, B and D genomes.",
                "Here we present a chromosome-scale assembly of T. aestivum cv. Chinese Spring.",
            ],
            "evidence_type": "explicit_statement",
            "stated_explicitly": True,
            "confidence": "high",
            "ambiguity_note": "Passage also mentions diploid (2x) progenitors; these are relatives, not the target.",
            "sex": "unknown",
            "strain_cultivar": "Chinese Spring",
        },
    },
    {
        "target": "Gossypium hirsutum",
        "passage": (
            "Upland cotton (Gossypium hirsutum) contains two distinct subgenomes, "
            "denoted A and D, derived from the merger of two diploid ancestors "
            "approximately 1-2 million years ago. We sequenced G. hirsutum acc. TM-1."
        ),
        "output": {
            "reasoning": (
                "G. hirsutum has two subgenomes (A and D) from two diploid ancestors. "
                "Two merged diploid genomes give a tetraploid. No explicit ploidy word "
                "or formula is given, so this is inference. The sequenced accession is TM-1."
            ),
            "target_organism": "Gossypium hirsutum",
            "ploidy_level": 4,
            "ploidy_category": "tetraploid",
            "mechanism": "allopolyploid",
            "chromosome_formula": None,
            "evidence_quotes": [
                "Upland cotton (Gossypium hirsutum) contains two distinct subgenomes, denoted A and D, derived from the merger of two diploid ancestors approximately 1-2 million years ago.",
                "We sequenced G. hirsutum acc. TM-1.",
            ],
            "evidence_type": "subgenome_count",
            "stated_explicitly": False,
            "confidence": "medium",
            "ambiguity_note": None,
            "sex": "unknown",
            "strain_cultivar": "TM-1",
        },
    },
    {
        "target": "Brassica cretica",
        "passage": (
            "We assembled the genome of Brassica cretica. The 250 Mb assembly was "
            "generated from PacBio HiFi reads and annotated with 28,000 "
            "protein-coding genes. Synteny with related crucifers was analyzed to "
            "study diversification."
        ),
        "output": {
            "reasoning": (
                "The passage reports assembly size, sequencing technology, gene count, "
                "and synteny but contains no ploidy-relevant evidence and no sex or "
                "strain information. I must not assume a ploidy from prior knowledge."
            ),
            "target_organism": "Brassica cretica",
            "ploidy_level": None,
            "ploidy_category": None,
            "mechanism": None,
            "chromosome_formula": None,
            "evidence_quotes": [],
            "evidence_type": "none",
            "stated_explicitly": False,
            "confidence": "low",
            "ambiguity_note": None,
            "sex": "unknown",
            "strain_cultivar": None,
        },
    },
]


def _user_block(target: str, passages_str: str, accession: str | None = None) -> str:
    target_clean = (
        target or ""
    ).strip() or "(not specified — use the organism whose genome this paper presents)"
    acc_line = f"Assembly accession: {accession}\n" if accession else ""
    return (
        f"Target organism: {target_clean}\n"
        f"{acc_line}\n"
        f"Passages:\n<passages>\n{passages_str.strip()}\n</passages>\n\n"
        "Return the single JSON object."
    )


def _build_user_prompt(species: str, passages: list, accession: str | None = None) -> str:
    blocks = []
    total = 0
    for i, passage in enumerate(passages, 1):
        text = (passage.get("text") or "").strip()
        if not text:
            continue
        if total + len(text) > MAX_CHARS:
            text = text[: max(0, MAX_CHARS - total)]
        blocks.append(f"[Passage {i} | section: {passage.get('section', '?')}]\n{text}")
        total += len(text)
        if total >= MAX_CHARS:
            break
    joined = "\n\n".join(blocks) if blocks else "(no passages)"
    return _user_block(species, joined, accession=accession)


def _build_messages(species: str, passages: list, accession: str | None = None) -> list:
    msgs = [{"role": "system", "content": _SYSTEM}]
    for ex in _FEWSHOT:
        msgs.append({"role": "user", "content": _user_block(str(ex["target"]), str(ex["passage"]))})
        msgs.append({"role": "assistant", "content": json.dumps(ex["output"], ensure_ascii=False)})
    msgs.append({"role": "user", "content": _build_user_prompt(species, passages, accession=accession)})
    return msgs


def _call_vllm(species: str, passages: list, accession: str | None = None) -> dict | None:
    payload = {
        "model": MODEL,
        "messages": _build_messages(species, passages, accession=accession),
        "temperature": 0,
        "max_tokens": 600,
        # vLLM extension: constrain output to the schema (ignored by plain OpenAI).
        "guided_json": _PLOIDY_SCHEMA,
        # llama-cpp-python / standard OpenAI fallback for JSON mode.
        "response_format": {"type": "json_object"},
    }
    try:
        r = requests.post(f"{BASE_URL}/chat/completions", json=payload, timeout=TIMEOUT)
        r.raise_for_status()
        content = r.json()["choices"][0]["message"]["content"]
    except (requests.RequestException, KeyError, ValueError, IndexError) as e:
        logger.warning(f"  [gemma] call failed ({e}); skipping Gemma layer")
        return None
    return _parse_json(content)


def _parse_json(text: str) -> dict | None:
    if not text:
        return None
    cleaned = text.strip()
    cleaned = re.sub(r"^```(?:json)?|```$", "", cleaned, flags=re.MULTILINE).strip()
    # grab the first {...} block if the model added prose
    match = re.search(r"\{.*\}", cleaned, flags=re.DOTALL)
    if match:
        cleaned = match.group()
    try:
        return json.loads(cleaned)
    except json.JSONDecodeError:
        return None


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", (text or "").lower()).strip()


def _is_grounded(quote: str, passages: list) -> bool:
    quote_norm = _norm(quote)
    if len(quote_norm) < 12:
        return False
    corpus = _norm(" ".join(passage.get("text", "") for passage in passages))
    if quote_norm in corpus:
        return True
    # Fuzzy fallback: ≥65% of content words (len>5) must appear in the corpus.
    # Using len>5 avoids short function words ("with", "from") skewing the score.
    words = [word for word in quote_norm.split() if len(word) > 5]
    if not words:
        return False
    return sum(1 for word in words if word in corpus) / len(words) >= 0.65


def _normalize_sex(raw_sex: object) -> str:
    """Map the model's free-form sex value onto the controlled vocabulary
    (_SEX_VALID / _SEX_MAP), defaulting to "unknown"."""
    sex_raw = str(raw_sex or "").lower().strip()
    return _SEX_MAP.get(sex_raw, sex_raw if sex_raw in _SEX_VALID else "unknown")


def _normalize_level(level: object) -> int | None:
    """Coerce the model's ploidy_level (int, numeric string, or ploidy word) to an
    integer. Rejects implausibly large values (> 20) — the model sometimes returns
    the 2n chromosome count (e.g. 30 from "2n=4x=30") instead of the ploidy."""
    if isinstance(level, str):
        level_word = level.lower().strip()
        if level_word in _WORD_LEVEL:
            level = _WORD_LEVEL[level_word]
        else:
            match = re.search(r"\d+", level_word)
            level = int(match.group()) if match else None
    if not isinstance(level, int) or level > 20:
        return None
    return level


def _normalize_mechanism(mechanism: object, quotes: list) -> str | None:
    """Map the model's mechanism onto the controlled VALID_MECHANISMS vocabulary,
    scanning for keywords when it returns a full ploidy word or a sentence, and
    inferring from the evidence quotes when it left the field null. Returns None
    if no controlled mechanism can be determined."""
    # Guard: model occasionally returns mechanism as a list (e.g. ["allopolyploid"])
    if isinstance(mechanism, list):
        mechanism = mechanism[0] if mechanism else None
    if isinstance(mechanism, str):
        mech_lower = mechanism.lower().strip()
        for keyword, controlled in _MECH_KEYWORDS:
            if keyword in mech_lower:
                mechanism = controlled
                break
        else:
            mechanism = None
    if not isinstance(mechanism, str) or mechanism not in VALID_MECHANISMS:
        mechanism = None

    # Infer mechanism from quotes when the model left it null
    if mechanism is None:
        all_quotes_text = " ".join(quotes).lower()
        for keyword, controlled in _MECH_KEYWORDS:
            if keyword in all_quotes_text:
                return controlled
    return mechanism


def _validate(obj: dict, passages: list) -> dict | None:
    if not isinstance(obj, dict):
        return None

    level = obj.get("ploidy_level")
    mechanism = obj.get("mechanism")
    explicit = bool(obj.get("stated_explicitly"))
    conf_str = str(obj.get("confidence") or "").lower().strip()
    chrom_formula = obj.get("chromosome_formula")
    evidence_type = obj.get("evidence_type")
    ambiguity = obj.get("ambiguity_note")

    # Handle evidence_quotes (list) with backward compat for legacy evidence_quote (str)
    raw_quotes = obj.get("evidence_quotes") or obj.get("evidence_quote") or []
    if isinstance(raw_quotes, str):
        raw_quotes = [raw_quotes] if raw_quotes else []
    quotes = [q.strip() for q in raw_quotes if q and isinstance(q, str)]
    quote = quotes[0] if quotes else ""  # primary quote for grounding / mechanism scan

    # Normalise the model's free-form fields onto controlled vocabularies.
    sex = _normalize_sex(obj.get("sex"))
    level = _normalize_level(level)
    mechanism = _normalize_mechanism(mechanism, quotes)

    # Strain / cultivar
    strain_cultivar = obj.get("strain_cultivar")
    if isinstance(strain_cultivar, str):
        strain_cultivar = strain_cultivar.strip() or None
    else:
        strain_cultivar = None

    # Grounding gate: at least one evidence_quote must be traceable to passages.
    grounded = any(_is_grounded(q, passages) for q in quotes) if quotes else False
    if (level is not None or mechanism is not None) and not grounded:
        logger.info("  [gemma] value rejected (evidence_quotes not grounded in passages)")
        level, mechanism = None, None

    # schema normalisation (mirror resolve_ploidy_fields)
    if level == 1:
        mechanism = None
    elif level is not None and level < 3:
        mechanism = None
    if mechanism == "amphidiploid" and level is not None and level % 2 != 0:
        mechanism = None  # amphidiploid implies an even level

    # Return None only if nothing useful: no ploidy AND no sex AND no strain
    has_sex = sex and sex not in ("unknown", "not_applicable")
    has_strain = strain_cultivar is not None
    if level is None and mechanism is None and not has_sex and not has_strain:
        return None

    # Map model's confidence string to a float (_CONF_MAP at module level).
    # Inferred results get an 88% multiplier — they're valuable but less certain.
    base_conf = _CONF_MAP.get(conf_str, 0.55)
    suggested = not explicit
    confidence = round(base_conf * (1.0 if explicit else 0.88), 2)
    source = "gemma_suggested" if suggested else "gemma_confirmed"
    method = "gemma_suggested" if suggested else "gemma"

    if suggested and level is not None:
        logger.info(f"  [gemma] suggested (inferred): level={level} mech={mechanism} conf={conf_str}")

    return {
        "level": level,
        "mechanism": mechanism,
        "irregular": False,
        "derivation": None,
        "tuple": [level, mechanism, False, None],
        "level_label": _LEVEL_LABEL.get(level) if level is not None else None,
        "status": "found_in_paper",
        "source": source,
        "confidence": confidence,
        "method": method,
        "suggested": suggested,
        "evidence": quote.strip()[:300],
        "evidence_quotes": [q[:300] for q in quotes],
        "evidence_type": evidence_type,
        "chromosome_formula": chrom_formula,
        "ambiguity_note": ambiguity,
        "sex": sex,
        "strain_cultivar": strain_cultivar,
    }


def run_gemma_ploidy(
    species: str, evidence_passages: list[dict], accession: str | None = None
) -> dict | None:
    """Run grounded ploidy/sex/strain extraction for `species` over the given
    evidence passages, returning a validated result dict (or None if the Gemma
    layer is disabled, there are no passages, or every backend fails). Backends
    are tried in order: HTTP server, in-process GGUF, then HF float32."""
    if not is_enabled():
        return None
    if not evidence_passages:
        return None
    logger.info("\n  [gemma] querying local Gemma for ploidy (grounded extraction) ...")
    # Priority: HTTP server > GGUF in-process > HF float32 direct
    if _server_is_up():
        raw = _call_vllm(species, evidence_passages, accession=accession)
        if raw is not None:
            return _validate(raw, evidence_passages)
    if _gguf_available():
        raw = _call_gguf(species, evidence_passages, accession=accession)
        if raw is not None:
            return _validate(raw, evidence_passages)
    if _weights_available():
        raw = _call_direct(species, evidence_passages, accession=accession)
        if raw is not None:
            return _validate(raw, evidence_passages)
    return None


def combine_with_gemma(ensemble: dict, gemma: dict | None) -> dict:
    """Fold a validated Gemma result into the rule+vector ensemble ploidy as a
    third vote. Inferred (`suggested`) results are attached as a non-destructive
    `gemma_suggestion`; explicit results can raise confidence on agreement, fill a
    missing level, or flag a conflict. Also backfills sex and cultivar when absent."""
    if not gemma:
        return ensemble

    out = dict(ensemble)
    ensemble_level = ensemble.get("level")
    gemma_level = gemma.get("level")
    suggested = gemma.get("suggested", False)

    # Annotation that is always attached when Gemma has something to say
    gemma_note = {
        "level": gemma_level,
        "mechanism": gemma.get("mechanism"),
        "evidence": gemma.get("evidence"),
        "suggested": suggested,
    }

    if suggested:
        # ── Gemma inferred (not explicitly stated) ──────────────────────────
        # Never overwrite the ensemble result; attach as a suggestion only.
        out["gemma_suggestion"] = gemma_note
        if ensemble_level is None and gemma_level is not None:
            # Nothing from rule/vector — show Gemma's inference as a low-confidence hint
            out["gemma_suggestion"]["note"] = (
                "No result from rule-based or vector search. "
                "Gemma inferred this level from genomic context — treat as a hint, "
                "not a confirmed answer."
            )
        elif ensemble_level is not None and gemma_level is not None and ensemble_level == gemma_level:
            # Gemma agrees with ensemble via inference — small confidence nudge
            out["confidence"] = min(round(ensemble.get("confidence", 0.5) + 0.05, 2), 0.95)
            out["gemma_suggestion"]["note"] = "Gemma inferred the same level; adds soft support."
        elif ensemble_level is not None and gemma_level is not None and ensemble_level != gemma_level:
            out["gemma_suggestion"]["note"] = (
                f"Gemma inferred level={gemma_level}, which conflicts with the ensemble "
                f"result of level={ensemble_level}. Manual review recommended."
            )
    else:
        # ── Gemma explicitly confirmed ───────────────────────────────────────
        if ensemble_level is not None and gemma_level is not None and ensemble_level == gemma_level:
            # All three sources agree → highest confidence
            out["confidence"] = min(round(ensemble.get("confidence", 0.5) + 0.2, 2), 0.98)
            out["method"] = "consensus_3way"
            if not out.get("mechanism") and gemma.get("mechanism"):
                out["mechanism"] = gemma["mechanism"]
                out["tuple"] = [
                    out["level"],
                    out["mechanism"],
                    out.get("irregular", False),
                    out.get("derivation"),
                ]
            out["gemma"] = gemma_note

        elif ensemble_level is None and gemma_level is not None:
            # Rule/vector found nothing; Gemma found an explicit statement
            out = dict(gemma)
            out["method"] = "gemma_only"
            out["confidence"] = round(gemma.get("confidence", 0.7) * 0.9, 2)

        elif ensemble_level is not None and gemma_level is not None and ensemble_level != gemma_level:
            # Explicit conflict — flag for manual review
            out["method"] = f"conflict_rule_vector_vs_gemma ({ensemble_level} vs {gemma_level})"
            out["confidence"] = round(ensemble.get("confidence", 0.5) * 0.8, 2)
            out["gemma"] = gemma_note

    # ── Apply sex from Gemma if ensemble has unknown ────────────────────────
    g_sex = gemma.get("sex")
    if g_sex and g_sex not in ("unknown", "not_applicable"):
        if out.get("sex", {}).get("value") in (None, "unknown"):
            out["sex"] = {"value": g_sex, "method": "gemma", "confidence": 0.65}

    # ── Add Gemma cultivar/strain if not already in the list ────────────────
    g_cultivar = gemma.get("strain_cultivar")
    if g_cultivar:
        cultivars = out.get("cultivars") or {}
        confirmed = cultivars.get("confirmed", []) or []
        unconfirmed = cultivars.get("unconfirmed", []) or []
        if g_cultivar not in confirmed and g_cultivar not in unconfirmed:
            out.setdefault("cultivars", {}).setdefault("unconfirmed", []).append(g_cultivar)

    return out
