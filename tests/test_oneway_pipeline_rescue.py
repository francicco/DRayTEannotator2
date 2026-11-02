"""Test curated-reference label transfer and safe policy changes end to end."""

import json
from pathlib import Path
import sys

from Bio import SeqIO
import pytest

from drayte.pipeline import classification as stage
from test_assisted_pipeline_rescue import assisted_case, legacy_outputs, read_table


@pytest.fixture
def many_query_case(assisted_case):
    case = assisted_case
    # Exercise the new external-reference default, independent of the older
    # reciprocal fixture used to test its backwards-compatible strict policy.
    case.config.extra["classification"].pop("mmseqs_rescue_mode", None)
    binary = Path(case.config.extra["classification"]["mmseqs_bin"])
    calls_path = case.reference.parent / "mmseqs_calls.jsonl"
    binary.write_text(
        f"#!{sys.executable}\n"
        "import json, sys\n"
        "from pathlib import Path\n"
        "def headers(path):\n"
        "    return [line[1:].split()[0] for line in Path(path).read_text().splitlines() if line.startswith('>')]\n"
        "query, target = headers(sys.argv[2]), headers(sys.argv[3])\n"
        "rows = []\n"
        "for item in query:\n"
        "    rows.append(f'{item}\\t{target[0]}\\t0.95\\t900\\t0\\t0\\t1\\t900\\t1\\t900\\t1e-80\\t1800\\n')\n"
        "Path(sys.argv[4]).write_text(''.join(rows))\n"
        f"with Path({str(calls_path)!r}).open('a') as handle:\n"
        "    handle.write(json.dumps({'query': query, 'target': target}) + '\\n')\n"
    )
    return case


def cache_state(case):
    return {path: (path.read_bytes(), path.stat().st_mtime_ns)
            for path in case.cached_paths}


def test_default_oneway_classifies_many_queries_from_one_reference(many_query_case):
    case = many_query_case
    original_sequences = {record.id.split("#")[0]: str(record.seq)
                          for record in SeqIO.parse(case.structure, "fasta")}
    before = cache_state(case)
    result = case.run()

    assert len(case.calls()) == 1
    assert case.calls()[0] == {
        "query": ["T-family-1#Unknown", "T-family-2#Unknown"],
        "target": ["Curated-Helitron#RC/Helitron"],
    }
    classifications = read_table(result["classifications_tsv"])
    assert set(classifications) == set(original_sequences)
    assert {row["order"] for row in classifications.values()} == {"Helitron"}
    evidence = read_table(result["evidence_tsv"])
    assert {row["rescue_target"] for row in evidence.values()} == {
        "Curated-Helitron#RC/Helitron"
    }
    output = list(SeqIO.parse(result["classified_library"], "fasta"))
    assert {record.id.split("#")[0]: str(record.seq) for record in output} == original_sequences
    assert {record.id for record in output} == {
        "T-family-1#RC/Helitron", "T-family-2#RC/Helitron"
    }
    assert result["n_classified"] == 2 and result["n_unknown"] == 0
    report = result["mmseqs_reference_rescue"]
    assert report["mode"] == "one-way"
    assert report["forward_hits"] == 2
    assert report["reverse_hits"] == report["reciprocal_pairs"] == 0
    assert report["families_with_rescue"] == 2
    assert report["reverse_tsv"] is None
    assert "mmseqs_reciprocal_rescue" not in result
    assert result["mmseqs_rescue_impact"]["previously_unknown_rescued"] == 2
    assert cache_state(case) == before


def test_completed_oneway_stage_is_reused_after_worker_changes(many_query_case, monkeypatch):
    case = many_query_case
    first = case.run()
    case.config.threads = 8
    case.config.extra["classification"].update(cpu=4, jobs=3, chunks=12)

    def unexpected(*args, **kwargs):
        pytest.fail("A complete matching one-way stage should skip evidence work")

    monkeypatch.setattr(stage, "run_reference_mmseqs_rescue", unexpected)
    monkeypatch.setattr(stage, "detect_complexity_from_fasta", unexpected)
    second = case.run()
    assert second["classified_library"] == first["classified_library"]
    assert len(case.calls()) == 1


def test_explicit_default_policy_reuses_implicit_default_stage(many_query_case, monkeypatch):
    case = many_query_case
    first = case.run()
    case.config.extra["classification"].update(
        mmseqs_rescue_mode="one-way", mmseqs_rescue_score_fraction=0.95
    )

    def unexpected(*args, **kwargs):
        pytest.fail("Explicit defaults should match the already completed implicit policy")

    monkeypatch.setattr(stage, "run_reference_mmseqs_rescue", unexpected)
    monkeypatch.setattr(stage, "detect_complexity_from_fasta", unexpected)
    second = case.run()
    assert second["classified_library"] == first["classified_library"]
    assert len(case.calls()) == 1


def test_switching_reciprocal_to_oneway_refreshes_labels_and_keeps_other_evidence(many_query_case):
    case = many_query_case
    before = cache_state(case)
    case.config.extra["classification"]["mmseqs_rescue_mode"] = "reciprocal"
    strict = case.run()
    assert strict["n_classified"] == 1 and strict["n_unknown"] == 1
    assert len(case.calls()) == 2
    assert strict["mmseqs_reciprocal_rescue"]["reciprocal_pairs"] == 1

    case.config.extra["classification"]["mmseqs_rescue_mode"] = "one-way"
    transferred = case.run()
    assert transferred["n_classified"] == 2 and transferred["n_unknown"] == 0
    assert len(case.calls()) == 3
    assert transferred["mmseqs_reference_rescue"]["mode"] == "one-way"
    assert "mmseqs_reciprocal_rescue" not in transferred
    assert cache_state(case) == before


def test_changed_competing_score_fraction_refreshes_reference_rescue(many_query_case):
    case = many_query_case
    case.run()
    case.config.extra["classification"]["mmseqs_rescue_score_fraction"] = 0.90
    result = case.run()
    assert len(case.calls()) == 2
    assert result["mmseqs_reference_rescue"]["filter_parameters"]["min_score_fraction"] == 0.90
    assert result["n_classified"] == 2


def test_version_one_reciprocal_stage_is_refreshed_to_default_oneway(many_query_case):
    case = many_query_case
    case.config.extra["classification"]["mmseqs_rescue_mode"] = "reciprocal"
    result = case.run()
    manifest = case.outdir / "classification.manifest.json"
    historical = json.loads(manifest.read_text())
    historical["input_signature"]["mmseqs_rescue_algorithm_version"] = 1
    historical["input_signature"].pop("mmseqs_rescue_mode", None)
    historical["input_signature"].pop("mmseqs_rescue_score_fraction", None)
    historical.pop("mmseqs_reference_rescue", None)
    historical["mmseqs_reciprocal_rescue"]["algorithm_version"] = 1
    historical["mmseqs_reciprocal_rescue"].pop("mode", None)
    manifest.write_text(json.dumps(historical))
    before = cache_state(case)
    case.config.extra["classification"].pop("mmseqs_rescue_mode")

    refreshed = case.run()
    assert refreshed["n_classified"] == 2 and result["n_classified"] == 1
    assert refreshed["mmseqs_reference_rescue"]["algorithm_version"] == stage.EXTERNAL_RESCUE_ALGORITHM_VERSION
    assert refreshed["mmseqs_reference_rescue"]["mode"] == "one-way"
    assert len(case.calls()) == 3
    assert cache_state(case) == before


def test_reference_content_changes_still_require_explicit_input_resolution(many_query_case):
    case = many_query_case
    case.run()
    case.reference.write_text(case.reference.read_text().replace("ACGTC", "TCGTC", 1))
    with pytest.raises(RuntimeError, match="classification.mmseqs_rescue_db contents"):
        case.run()
    assert len(case.calls()) == 1


def test_rescue_policy_refresh_does_not_approve_other_biological_setting_changes(many_query_case):
    case = many_query_case
    case.run()
    case.config.extra["classification"].update(
        mmseqs_rescue_mode="reciprocal", min_orf_nt=600
    )
    with pytest.raises(RuntimeError, match="classification.min_orf_nt"):
        case.run()
    assert len(case.calls()) == 1


def test_unassisted_legacy_stage_does_not_require_external_policy(many_query_case, monkeypatch):
    case = many_query_case
    case.config.extra["classification"].pop("mmseqs_rescue_db")
    legacy_outputs(case)

    def unexpected(*args, **kwargs):
        pytest.fail("Unassisted matching completed stages should remain reusable")

    monkeypatch.setattr(stage, "run_reference_mmseqs_rescue", unexpected)
    monkeypatch.setattr(stage, "detect_complexity_from_fasta", unexpected)
    result = case.run()
    assert result["n_unknown"] == 2
    assert case.calls() == []
