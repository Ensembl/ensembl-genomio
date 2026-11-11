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
"""Unit testing of `ensembl.io.genomio.literature.gemma` module.

The language model is never loaded: `llama_cpp`, `torch` and `transformers` are replaced with mocks.

Typical usage example::
    $ pytest test_gemma.py

"""

# pylint: disable=protected-access,too-many-positional-arguments

import json
import os
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
import requests

from ensembl.io.genomio.literature import gemma

PASSAGES = [{"section": "Abstract", "text": "Peanut (Arachis hypogaea) is an allotetraploid, 2n = 4x = 40."}]
QUOTE = "Peanut (Arachis hypogaea) is an allotetraploid"


@pytest.mark.parametrize(
    ("file_name", "create", "expected"),
    [
        pytest.param("model.gguf", True, True, id="GGUF file present"),
        pytest.param("model.gguf", False, False, id="GGUF file missing"),
        pytest.param("model.bin", True, False, id="Not a GGUF file"),
    ],
)
def test_gguf_available(tmp_path: Path, file_name: str, create: bool, expected: bool) -> None:
    """Test :func:`gemma._gguf_available`.

    Args:
        tmp_path: Pytest-provided temporary directory.
        file_name: Configured model file name.
        create: Whether the file exists.
        expected: Expected availability.

    """
    model = tmp_path / file_name
    if create:
        model.touch()
    with patch.object(gemma, "MODEL", str(model)):
        assert gemma._gguf_available() == expected


def test_weights_available(tmp_path: Path) -> None:
    """Test :func:`gemma._weights_available`.

    Args:
        tmp_path: Pytest-provided temporary directory.

    """
    with patch.object(gemma, "_HF_MODEL_DIR", str(tmp_path)):
        assert not gemma._weights_available()
        (tmp_path / "config.json").touch()
        assert gemma._weights_available()


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param(None, True, id="Server answers"),
        pytest.param(requests.ConnectionError("down"), False, id="No server"),
    ],
)
@patch("ensembl.io.genomio.literature.gemma.requests.get")
def test_server_is_up(mock_get: MagicMock, side_effect: Exception | None, expected: bool) -> None:
    """Test :func:`gemma._server_is_up`.

    Args:
        mock_get: Patched `requests.get()`.
        side_effect: Exception raised by the request, if any.
        expected: Expected server state.

    """
    mock_get.side_effect = side_effect
    assert gemma._server_is_up() == expected


@pytest.mark.parametrize(
    ("env", "detected", "expected"),
    [
        pytest.param({"GEMMA_ENABLED": "1"}, (False, False, False), True, id="Forced on"),
        pytest.param({"GEMMA_ENABLED": "0"}, (True, True, True), False, id="Forced off"),
        pytest.param({}, (False, False, True), True, id="Server detected"),
        pytest.param({}, (False, False, False), False, id="Nothing detected"),
    ],
)
@patch.dict("os.environ")
@patch("ensembl.io.genomio.literature.gemma._server_is_up")
@patch("ensembl.io.genomio.literature.gemma._weights_available")
@patch("ensembl.io.genomio.literature.gemma._gguf_available")
def test_is_enabled(
    mock_gguf: MagicMock,
    mock_weights: MagicMock,
    mock_server: MagicMock,
    env: dict[str, str],
    detected: tuple[bool, bool, bool],
    expected: bool,
) -> None:
    """Test :func:`gemma.is_enabled`.

    Args:
        mock_gguf: Patched `_gguf_available()`.
        mock_weights: Patched `_weights_available()`.
        mock_server: Patched `_server_is_up()`.
        env: ``GEMMA_ENABLED`` setting, if any.
        detected: What the GGUF, weights and server checks report.
        expected: Expected state.

    """
    os.environ.pop("GEMMA_ENABLED", None)
    os.environ.update(env)
    mock_gguf.return_value, mock_weights.return_value, mock_server.return_value = detected
    assert gemma.is_enabled() == expected


@patch("ensembl.io.genomio.literature.gemma._gguf_model", None)
def test_load_gguf_backend() -> None:
    """Test :func:`gemma._load_gguf_backend` loads the model only once."""
    llama_cpp = MagicMock()
    with patch.dict("sys.modules", {"llama_cpp": llama_cpp}):
        gemma._load_gguf_backend()
        gemma._load_gguf_backend()
    assert gemma._gguf_model is llama_cpp.Llama.return_value
    llama_cpp.Llama.assert_called_once()


ANSWER = {"choices": [{"message": {"content": '{"ploidy_level": 4}'}}]}


@pytest.mark.parametrize(
    ("load_error", "model", "expected"),
    [
        pytest.param(None, MagicMock(**{"create_chat_completion.return_value": ANSWER}), {"ploidy_level": 4}),
        pytest.param(ImportError("no llama_cpp"), MagicMock(), None, id="Load fails"),
        pytest.param(None, None, None, id="Model not initialised"),
        pytest.param(
            None, MagicMock(**{"create_chat_completion.side_effect": RuntimeError}), None, id="Fails"
        ),
    ],
)
@patch("ensembl.io.genomio.literature.gemma._load_gguf_backend")
def test_call_gguf(
    mock_load: MagicMock,
    load_error: Exception | None,
    model: MagicMock | None,
    expected: dict | None,
) -> None:
    """Test :func:`gemma._call_gguf`.

    Args:
        mock_load: Patched `_load_gguf_backend()`.
        load_error: Exception raised while loading, if any.
        model: Loaded GGUF model.
        expected: Expected parsed answer.

    """
    mock_load.side_effect = load_error
    with patch.object(gemma, "_gguf_model", model):
        assert gemma._call_gguf("Arachis hypogaea", PASSAGES) == expected


@patch("ensembl.io.genomio.literature.gemma._direct_tokenizer", None)
@patch("ensembl.io.genomio.literature.gemma._direct_model", None)
def test_load_direct_backend() -> None:
    """Test :func:`gemma._load_direct_backend` loads the model and tokenizer only once."""
    transformers = MagicMock()
    with patch.dict("sys.modules", {"torch": MagicMock(), "transformers": transformers}):
        gemma._load_direct_backend()
        gemma._load_direct_backend()
    assert (
        gemma._direct_model
        is transformers.AutoModelForCausalLM.from_pretrained.return_value.eval.return_value
    )
    assert gemma._direct_tokenizer is transformers.AutoTokenizer.from_pretrained.return_value
    transformers.AutoModelForCausalLM.from_pretrained.assert_called_once()


INPUT_IDS = MagicMock(spec=["shape"], shape=(1, 3))
MODEL = MagicMock(**{"generate.return_value": [[5, 6, 7, 8, 9]]})


@pytest.mark.parametrize(
    ("load_error", "model", "encoded", "eot", "expected_eos", "expected"),
    [
        pytest.param(
            None,
            MODEL,
            MagicMock(input_ids=INPUT_IDS, attention_mask="mask"),
            106,
            [1, 106],
            {"ploidy_level": 4},
            id="Batch encoding with end-of-turn token",
        ),
        pytest.param(None, MODEL, INPUT_IDS, None, [1], {"ploidy_level": 4}, id="Plain tensor"),
        pytest.param(OSError("no weights"), MODEL, INPUT_IDS, None, None, None, id="Load fails"),
        pytest.param(None, None, INPUT_IDS, None, None, None, id="Model not initialised"),
        pytest.param(
            None,
            MagicMock(**{"generate.side_effect": RuntimeError("OOM")}),
            INPUT_IDS,
            None,
            None,
            None,
            id="Generate fails",
        ),
    ],
)
@patch("ensembl.io.genomio.literature.gemma._load_direct_backend")
def test_call_direct(
    mock_load: MagicMock,
    load_error: Exception | None,
    model: MagicMock | None,
    encoded: MagicMock,
    eot: int | None,
    expected_eos: list[int] | None,
    expected: dict | None,
) -> None:
    """Test :func:`gemma._call_direct` stops on the right tokens and decodes only the new tokens.

    Args:
        mock_load: Patched `_load_direct_backend()`.
        load_error: Exception raised while loading, if any.
        model: Loaded model.
        encoded: What the chat template returns.
        eot: ID of the ``<end_of_turn>`` token, or None when the tokenizer lacks it.
        expected_eos: Expected stop-token IDs, if generation runs.
        expected: Expected parsed answer.

    """
    mock_load.side_effect = load_error
    tokenizer = MagicMock(eos_token_id=1, **{"apply_chat_template.return_value": encoded})
    tokenizer.convert_tokens_to_ids.return_value = eot
    tokenizer.decode.return_value = '{"ploidy_level": 4}'
    with (
        patch.dict("sys.modules", {"torch": MagicMock()}),
        patch.object(gemma, "_direct_model", model),
        patch.object(gemma, "_direct_tokenizer", tokenizer),
    ):
        assert gemma._call_direct("Arachis hypogaea", PASSAGES) == expected
    if expected is not None:
        assert model.generate.call_args.kwargs["eos_token_id"] == expected_eos  # type: ignore[union-attr]
        tokenizer.decode.assert_called_once_with([8, 9], skip_special_tokens=True)


@pytest.mark.parametrize(
    ("side_effect", "expected"),
    [
        pytest.param([MagicMock(**{"json.return_value": ANSWER})], {"ploidy_level": 4}, id="Answer parsed"),
        pytest.param([MagicMock(**{"json.return_value": {"choices": []}})], None, id="No choices"),
        pytest.param(requests.Timeout("slow"), None, id="Request error"),
    ],
)
@patch("ensembl.io.genomio.literature.gemma.requests.post")
def test_call_vllm(
    mock_post: MagicMock,
    side_effect: Any,
    expected: dict | None,
) -> None:
    """Test :func:`gemma._call_vllm`.

    Args:
        mock_post: Patched `requests.post()`.
        side_effect: Response returned, or exception raised, by the server.
        expected: Expected parsed answer.

    """
    mock_post.side_effect = side_effect
    assert gemma._call_vllm("Arachis hypogaea", PASSAGES) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        pytest.param('{"ploidy_level": 4}', {"ploidy_level": 4}, id="Plain JSON"),
        pytest.param('```json\n{"ploidy_level": 4}\n```', {"ploidy_level": 4}, id="Code fence"),
        pytest.param(
            'Here you go: {"ploidy_level": 4} Hope it helps.', {"ploidy_level": 4}, id="Prose around"
        ),
        pytest.param("", None, id="Empty"),
        pytest.param("{ploidy: 4}", None, id="Invalid JSON"),
    ],
)
def test_parse_json(text: str, expected: dict | None) -> None:
    """Test :func:`gemma._parse_json`.

    Args:
        text: Raw model output.
        expected: Expected parsed object.

    """
    assert gemma._parse_json(text) == expected


@pytest.mark.parametrize(
    ("quote", "expected"),
    [
        pytest.param("4x = 40.", False, id="Too short"),
        pytest.param("PEANUT (Arachis   hypogaea) is an", True, id="Exact after normalising"),
        pytest.param("Arachis hypogaea, an allotetraploid peanut", True, id="Most content words present"),
        pytest.param("Arachis hypogaea shows strong drought tolerance", False, id="Too few content words"),
        pytest.param("it is a 4x one, as we see", False, id="No content words"),
    ],
)
def test_is_grounded(quote: str, expected: bool) -> None:
    """Test :func:`gemma._is_grounded`.

    Args:
        quote: Evidence quote returned by the model.
        expected: Whether the quote is grounded in the passages.

    """
    assert gemma._is_grounded(quote, PASSAGES) == expected


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        pytest.param(" Female ", "female", id="Controlled term"),
        pytest.param("hermaphroditic", "hermaphrodite", id="Mapped synonym"),
        pytest.param("boy", "unknown", id="Not in vocabulary"),
        pytest.param(None, "unknown", id="Missing"),
    ],
)
def test_normalize_sex(raw: object, expected: str) -> None:
    """Test :func:`gemma._normalize_sex`.

    Args:
        raw: Model's sex value.
        expected: Expected controlled term.

    """
    assert gemma._normalize_sex(raw) == expected


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        pytest.param(4, 4, id="Integer"),
        pytest.param(" Hexaploid", 6, id="Ploidy word"),
        pytest.param("8x", 8, id="Number in a string"),
        pytest.param(30, None, id="Chromosome count, not a level"),
        pytest.param("polyploid", None, id="No level"),
        pytest.param(None, None, id="Missing"),
    ],
)
def test_normalize_level(raw: object, expected: int | None) -> None:
    """Test :func:`gemma._normalize_level`.

    Args:
        raw: Model's ploidy level.
        expected: Expected integer level.

    """
    assert gemma._normalize_level(raw) == expected


@pytest.mark.parametrize(
    ("mechanism", "quotes", "expected"),
    [
        pytest.param("autopolyploid", [], "autopolyploid", id="Controlled term"),
        pytest.param(
            ["Segmental allopolyploid"], [], "segmental_allopolyploid", id="List, most specific keyword"
        ),
        pytest.param("an allotetraploid hybrid", [], "allopolyploid", id="Keyword in a sentence"),
        pytest.param("hybrid", ["It is an amphidiploid."], "amphidiploid", id="Inferred from quotes"),
        pytest.param(None, ["A diploid."], None, id="Nothing found"),
    ],
)
def test_normalize_mechanism(mechanism: object, quotes: list[str], expected: str | None) -> None:
    """Test :func:`gemma._normalize_mechanism`.

    Args:
        mechanism: Model's mechanism value.
        quotes: Evidence quotes.
        expected: Expected controlled mechanism.

    """
    assert gemma._normalize_mechanism(mechanism, quotes) == expected


@pytest.mark.parametrize(
    ("answer", "expected"),
    [
        pytest.param(
            {
                "ploidy_level": "tetraploid",
                "mechanism": "allotetraploid",
                "stated_explicitly": True,
                "confidence": "High",
                "evidence_quotes": [f"  {QUOTE} ", None],
                "sex": "unknown",
                "strain_cultivar": " Tifrunner ",
            },
            (4, "allopolyploid", 0.78, "gemma_confirmed", "unknown", "Tifrunner"),
            id="Explicit and grounded",
        ),
        pytest.param("not a dict", None, id="Not a dict"),
        pytest.param(
            {"ploidy_level": 4, "evidence_quotes": ["Unrelated sentence about rice"]}, None, id="Ungrounded"
        ),
        pytest.param(
            {"ploidy_level": 4, "evidence_quotes": ["Unrelated sentence about rice"], "sex": "female"},
            (None, None, 0.48, "gemma_suggested", "female", None),
            id="Ungrounded level dropped, sex kept",
        ),
        pytest.param(
            {"ploidy_level": 4, "confidence": "medium", "evidence_quote": QUOTE},
            (4, "allopolyploid", 0.55, "gemma_suggested", "unknown", None),
            id="Inferred, legacy quote field",
        ),
        pytest.param(
            {
                "ploidy_level": 2,
                "mechanism": "allopolyploid",
                "stated_explicitly": True,
                "evidence_quotes": [QUOTE],
            },
            (2, None, 0.55, "gemma_confirmed", "unknown", None),
            id="Diploid has no mechanism",
        ),
        pytest.param(
            {
                "ploidy_level": 3,
                "mechanism": "amphidiploid",
                "stated_explicitly": True,
                "evidence_quotes": [QUOTE],
            },
            (3, None, 0.55, "gemma_confirmed", "unknown", None),
            id="Amphidiploid needs an even level",
        ),
        pytest.param(
            {"ploidy_level": None, "strain_cultivar": "TM-1", "sex": "not_applicable"},
            (None, None, 0.48, "gemma_suggested", "not_applicable", "TM-1"),
            id="Strain only",
        ),
    ],
)
def test_validate(answer: Any, expected: tuple | None) -> None:
    """Test :func:`gemma._validate` grounding, schema normalisation and confidence.

    Args:
        answer: Parsed model answer.
        expected: Expected ``(level, mechanism, confidence, source, sex, strain_cultivar)``, or None.

    """
    result = gemma._validate(answer, PASSAGES)
    fields = ("level", "mechanism", "confidence", "source", "sex", "strain_cultivar")
    assert (tuple(result[field] for field in fields) if result else None) == expected


def test_user_block() -> None:
    """Test :func:`gemma._user_block`."""
    assert gemma._user_block(" Arachis hypogaea ", "  text  ", accession="GCA_1") == (
        "Target organism: Arachis hypogaea\nAssembly accession: GCA_1\n\n"
        "Passages:\n<passages>\ntext\n</passages>\n\nReturn the single JSON object."
    )
    assert gemma._user_block("", "text").startswith(
        "Target organism: (not specified — use the organism whose genome this paper presents)\n\nPassages:"
    )


@pytest.mark.parametrize(
    ("passages", "expected"),
    [
        pytest.param(
            [{"text": "abcd", "section": "Abstract"}, {"text": " "}, {"text": "efghij"}, {"text": "klm"}],
            "[Passage 1 | section: Abstract]\nabcd\n\n[Passage 3 | section: ?]\nefgh",
            id="Empty skipped, truncated at the cap",
        ),
        pytest.param([], "(no passages)", id="No passages"),
    ],
)
@patch("ensembl.io.genomio.literature.gemma.MAX_CHARS", 8)
@patch("ensembl.io.genomio.literature.gemma._user_block", side_effect=lambda target, text, accession: text)
def test_build_user_prompt(mock_block: MagicMock, passages: list[dict], expected: str) -> None:
    """Test :func:`gemma._build_user_prompt`.

    Args:
        mock_block: Patched `_user_block()`, returning the passage text it receives.
        passages: Evidence passages.
        expected: Expected passage block.

    """
    assert gemma._build_user_prompt("Arachis hypogaea", passages, accession="GCA_1") == expected
    mock_block.assert_called_once_with("Arachis hypogaea", expected, accession="GCA_1")


@patch("ensembl.io.genomio.literature.gemma._build_user_prompt", return_value="final prompt")
def test_build_messages(mock_prompt: MagicMock) -> None:
    """Test :func:`gemma._build_messages`.

    Args:
        mock_prompt: Patched `_build_user_prompt()`.

    """
    messages = gemma._build_messages("Arachis hypogaea", PASSAGES, accession="GCA_1")
    assert [message["role"] for message in messages] == ["system"] + ["user", "assistant"] * 3 + ["user"]
    assert messages[0]["content"] == gemma._SYSTEM
    assert messages[1]["content"] == gemma._user_block("Triticum aestivum", str(gemma._FEWSHOT[0]["passage"]))
    assert json.loads(messages[2]["content"]) == gemma._FEWSHOT[0]["output"]
    assert messages[-1]["content"] == "final prompt"
    mock_prompt.assert_called_once_with("Arachis hypogaea", PASSAGES, accession="GCA_1")


@pytest.mark.parametrize(
    ("enabled", "passages", "available", "answers", "expected_backend"),
    [
        pytest.param(False, PASSAGES, (True, True, True), ("vllm", "gguf", "direct"), None, id="Disabled"),
        pytest.param(True, [], (True, True, True), ("vllm", "gguf", "direct"), None, id="No passages"),
        pytest.param(
            True, PASSAGES, (True, True, True), ("vllm", "gguf", "direct"), "vllm", id="Server first"
        ),
        pytest.param(True, PASSAGES, (True, True, True), (None, "gguf", "direct"), "gguf", id="Server fails"),
        pytest.param(
            True, PASSAGES, (False, False, True), (None, None, "direct"), "direct", id="Weights only"
        ),
        pytest.param(True, PASSAGES, (True, True, True), (None, None, None), None, id="Every backend fails"),
    ],
)
@patch("ensembl.io.genomio.literature.gemma._validate", lambda raw, passages: {"from": raw})
@patch("ensembl.io.genomio.literature.gemma._call_direct")
@patch("ensembl.io.genomio.literature.gemma._call_gguf")
@patch("ensembl.io.genomio.literature.gemma._call_vllm")
@patch("ensembl.io.genomio.literature.gemma._weights_available")
@patch("ensembl.io.genomio.literature.gemma._gguf_available")
@patch("ensembl.io.genomio.literature.gemma._server_is_up")
@patch("ensembl.io.genomio.literature.gemma.is_enabled")
def test_run_gemma_ploidy(
    mock_enabled: MagicMock,
    mock_server: MagicMock,
    mock_gguf: MagicMock,
    mock_weights: MagicMock,
    mock_vllm: MagicMock,
    mock_call_gguf: MagicMock,
    mock_direct: MagicMock,
    enabled: bool,
    passages: list[dict],
    available: tuple[bool, bool, bool],
    answers: tuple,
    expected_backend: str | None,
) -> None:
    """Test :func:`gemma.run_gemma_ploidy` tries the backends in order and validates the first answer.

    Args:
        mock_enabled: Patched `is_enabled()`.
        mock_server: Patched `_server_is_up()`.
        mock_gguf: Patched `_gguf_available()`.
        mock_weights: Patched `_weights_available()`.
        mock_vllm: Patched `_call_vllm()`.
        mock_call_gguf: Patched `_call_gguf()`.
        mock_direct: Patched `_call_direct()`.
        enabled: Whether the Gemma layer is enabled.
        passages: Evidence passages.
        available: Whether the server, GGUF file and HF weights are available.
        answers: Answer returned by the server, GGUF and HF backends.
        expected_backend: Backend whose answer is expected to be returned.

    """
    mock_enabled.return_value = enabled
    mock_server.return_value, mock_gguf.return_value, mock_weights.return_value = available
    mock_vllm.return_value, mock_call_gguf.return_value, mock_direct.return_value = answers
    result = gemma.run_gemma_ploidy("Arachis hypogaea", passages)
    assert result == ({"from": expected_backend} if expected_backend else None)


ENSEMBLE = {
    "level": 4,
    "mechanism": None,
    "irregular": False,
    "derivation": None,
    "tuple": [4, None, False, None],
    "confidence": 0.7,
    "method": "consensus",
}


@pytest.mark.parametrize(
    ("ensemble", "answer", "expected"),
    [
        pytest.param(ENSEMBLE, None, ENSEMBLE, id="No Gemma result"),
        pytest.param(
            ENSEMBLE,
            {"level": 6, "suggested": True, "sex": "unknown"},
            {"method": "consensus", "confidence": 0.7},
            id="Inferred answer discarded",
        ),
        pytest.param(
            ENSEMBLE,
            {"level": 4, "mechanism": "allopolyploid", "evidence": "e", "suggested": False},
            {
                "method": "consensus_3way",
                "confidence": 0.9,
                "mechanism": "allopolyploid",
                "tuple": [4, "allopolyploid", False, None],
                "gemma": {"level": 4, "mechanism": "allopolyploid", "evidence": "e", "suggested": False},
            },
            id="Agreement",
        ),
        pytest.param(
            {**ENSEMBLE, "confidence": 0.9, "mechanism": "autopolyploid"},
            {"level": 4, "mechanism": "allopolyploid"},
            {
                "method": "consensus_3way",
                "confidence": 0.98,
                "mechanism": "autopolyploid",
                "gemma": {"level": 4, "mechanism": "allopolyploid", "evidence": None, "suggested": False},
            },
            id="Agreement, capped, mechanism kept",
        ),
        pytest.param(
            {**ENSEMBLE, "level": None},
            {
                "level": 6,
                "confidence": 0.78,
                "source": "gemma_confirmed",
                "sex": "female",
                "strain_cultivar": "x",
            },
            {"level": 6, "method": "gemma_only", "confidence": 0.7, "source": "gemma_confirmed"},
            id="Gemma fills a missing level",
        ),
        pytest.param(
            ENSEMBLE,
            {"level": 2, "evidence": "e"},
            {
                "method": "conflict_rule_vector_vs_gemma (4 vs 2)",
                "confidence": 0.56,
                "gemma": {"level": 2, "mechanism": None, "evidence": "e", "suggested": False},
            },
            id="Conflict flagged",
        ),
    ],
)
def test_combine_with_gemma(ensemble: dict, answer: dict | None, expected: dict) -> None:
    """Test :func:`gemma.combine_with_gemma` combines the ploidy only, never Gemma's sex or cultivar.

    Args:
        ensemble: Rule + vector ensemble ploidy result.
        answer: Validated Gemma result, if any.
        expected: Fields expected in the combined result.

    """
    result = gemma.combine_with_gemma(dict(ensemble), answer)
    assert expected.items() <= result.items()
    assert "sex" not in result and "strain_cultivar" not in result


RECORD = {"sex": {"value": "unknown"}, "cultivars": {"confirmed": ["TM-1"], "unconfirmed": ["x"]}}


@pytest.mark.parametrize(
    ("record", "answer", "expected"),
    [
        pytest.param(RECORD, None, RECORD, id="No Gemma result"),
        pytest.param(
            RECORD,
            {"sex": "female", "strain_cultivar": "Tifrunner"},
            {
                "sex": {"value": "female", "method": "gemma", "confidence": 0.65},
                "cultivars": {"confirmed": ["TM-1"], "unconfirmed": ["x", "Tifrunner"]},
            },
            id="Sex and cultivar backfilled",
        ),
        pytest.param(
            {"sex": {"value": "male"}, "cultivars": {"confirmed": ["TM-1"]}},
            {"sex": "female", "strain_cultivar": "TM-1"},
            {"sex": {"value": "male"}, "cultivars": {"confirmed": ["TM-1"]}},
            id="Known sex and cultivar kept",
        ),
        pytest.param({}, {"sex": "unknown"}, {}, id="Gemma has no sex either"),
    ],
)
def test_backfill_from_gemma(record: dict, answer: dict | None, expected: dict) -> None:
    """Test :func:`gemma.backfill_from_gemma` fills only what the record is missing, without changing it.

    Args:
        record: Final pipeline record.
        answer: Validated Gemma result, if any.
        expected: Expected record.

    """
    before = json.dumps(record)
    assert gemma.backfill_from_gemma(record, answer) == expected
    assert json.dumps(record) == before
