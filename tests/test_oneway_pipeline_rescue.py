"""Test curated-reference label transfer and safe policy changes end to end."""

import json
import csv
import hashlib
from pathlib import Path
import sys

from Bio import SeqIO
import pytest

from drayte.pipeline import classification as stage
from test_assisted_pipeline_rescue import assisted_case, legacy_outputs, read_table


@pytest.fixture
def many_query_case(assisted_case):
    case = assisted_case
    # Exercise the external-reference default without legacy configuration.
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
    assert "reverse_hits" not in report and "reciprocal_pairs" not in report
    assert report["families_with_rescue"] == 2
    assert "reverse_tsv" not in report
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


def seed_historical_reciprocal_result(case, *, version, report_key):
    """Create a realistic historical result without executing removed code."""
    current = case.run()
    manifest = case.outdir / "classification.manifest.json"
    historical = json.loads(manifest.read_text())
    historical["n_classified"] = 1
    historical["n_unknown"] = 1
    historical["input_signature"]["mmseqs_rescue_algorithm_version"] = version
    if version == 1:
        historical["input_signature"].pop("mmseqs_rescue_mode", None)
        historical["input_signature"].pop("mmseqs_rescue_score_fraction", None)
    else:
        historical["input_signature"]["mmseqs_rescue_mode"] = "reciprocal"
    report = historical.pop("mmseqs_reference_rescue")
    report.update(algorithm_version=version, mode="reciprocal", reverse_hits=1,
                  reciprocal_pairs=1, families_with_rescue=1, rescued_families=1)
    historical[report_key] = report
    classified = Path(current["classified_library"])
    classified.write_text(classified.read_text().replace(
        ">T-family-2#RC/Helitron", ">T-family-2#Unknown"
    ))
    for key in ("classifications_tsv", "evidence_tsv"):
        path = Path(current[key])
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames
            rows = list(reader)
        for row in rows:
            if row["family_id"] == "T-family-2":
                for field in ("class", "order", "superfamily", "final_class",
                              "final_order", "final_superfamily", "rescue_class",
                              "rescue_order", "rescue_superfamily"):
                    if field in row:
                        row[field] = "Unknown"
                if "status" in row:
                    row["status"] = "unknown"
        with path.open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    manifest.write_text(json.dumps(historical))
    checkpoint = case.outdir / "mmseqs/mmseqs_rescue.checkpoint.json"
    completed = json.loads(checkpoint.read_text())
    completed["input_signature"].update(algorithm_version=version, mode="reciprocal")
    reverse = case.outdir / "mmseqs/mmseqs_rescue.reverse.tsv"
    reverse.write_text(
        "Curated-Helitron#RC/Helitron\tT-family-1#Unknown\t0.95\t900\t0\t0\t1\t900\t1\t900\t1e-80\t1800\n"
    )
    completed["reverse_sha256"] = hashlib.sha256(reverse.read_bytes()).hexdigest()
    checkpoint.write_text(json.dumps(completed))
    Path(current["mmseqs_rescue_report"]).write_text(json.dumps(report))
    (case.reference.parent / "mmseqs_calls.jsonl").write_text("")
    return historical


@pytest.mark.parametrize("version", [1, 2])
@pytest.mark.parametrize("report_key", ["mmseqs_reciprocal_rescue", "mmseqs_reference_rescue"])
def test_historical_reciprocal_result_refreshes_labels_and_keeps_other_evidence(
        many_query_case, version, report_key, caplog):
    case = many_query_case
    before = cache_state(case)
    historical = seed_historical_reciprocal_result(case, version=version, report_key=report_key)
    assert historical["n_classified"] == 1 and historical["n_unknown"] == 1
    assert case.calls() == []
    # Existing cluster YAML should migrate without requiring a fresh outdir.
    case.config.extra["classification"]["mmseqs_rescue_mode"] = "reciprocal"
    transferred = case.run()
    assert transferred["n_classified"] == 2 and transferred["n_unknown"] == 0
    assert len(case.calls()) == 1
    assert transferred["mmseqs_reference_rescue"]["mode"] == "one-way"
    assert "mmseqs_reciprocal_rescue" not in transferred
    assert "reverse_hits" not in transferred["mmseqs_reference_rescue"]
    assert "reciprocal_pairs" not in transferred["mmseqs_reference_rescue"]
    assert "Reciprocal rescue has been removed" in caplog.text
    assert cache_state(case) == before


def test_changed_competing_score_fraction_refreshes_reference_rescue(many_query_case):
    case = many_query_case
    case.run()
    case.config.extra["classification"]["mmseqs_rescue_score_fraction"] = 0.90
    result = case.run()
    assert len(case.calls()) == 2
    assert result["mmseqs_reference_rescue"]["filter_parameters"]["min_score_fraction"] == 0.90
    assert result["n_classified"] == 2


def test_historical_completed_oneway_result_remains_reusable(many_query_case, monkeypatch):
    case = many_query_case
    first = case.run()
    before = cache_state(case)
    # Already completed version-two one-way runs also carried null/zero fields
    # for the formerly available reverse search. Those outputs remain valid.
    manifest = case.outdir / "classification.manifest.json"
    historical = json.loads(manifest.read_text())
    historical["mmseqs_reference_rescue"].update(
        reverse_tsv=None, reverse_sha256=None, reverse_hits=0, reciprocal_pairs=0
    )
    manifest.write_text(json.dumps(historical))

    def unexpected(*args, **kwargs):
        pytest.fail("Removing the reverse implementation should preserve valid one-way results")

    monkeypatch.setattr(stage, "run_reference_mmseqs_rescue", unexpected)
    monkeypatch.setattr(stage, "detect_complexity_from_fasta", unexpected)
    second = case.run()
    assert second["classified_library"] == first["classified_library"]
    assert second["n_classified"] == 2
    assert len(case.calls()) == 1
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
    (case.outdir / "mmseqs/mmseqs_rescue.tsv").unlink()

    def unexpected(*args, **kwargs):
        pytest.fail("Unassisted matching completed stages should remain reusable")

    monkeypatch.setattr(stage, "run_reference_mmseqs_rescue", unexpected)
    monkeypatch.setattr(stage, "detect_complexity_from_fasta", unexpected)
    result = case.run()
    assert result["n_unknown"] == 2
    assert case.calls() == []


def test_legacy_within_library_tsv_refreshes_labels_once_without_external_search(
        many_query_case, monkeypatch):
    case = many_query_case
    case.config.extra["classification"].pop("mmseqs_rescue_db")
    legacy_outputs(case)
    legacy_tsv = case.outdir / "mmseqs/mmseqs_rescue.tsv"
    legacy_tsv.write_text(
        "T-family-1#Unknown\tT-family-2#RC/Helitron\t0.95\t900\t0\t0\t1\t900\t1\t900\t1e-80\t1800\n"
    )
    before = cache_state(case)

    def no_external_search(*args, **kwargs):
        pytest.fail("Within-library TSV evidence should not run an external search")

    monkeypatch.setattr(stage, "run_reference_mmseqs_rescue", no_external_search)
    first = case.run()
    classified = read_table(first["classifications_tsv"])
    assert classified["T-family-1"]["order"] == "Helitron"
    assert classified["T-family-2"]["order"] == "Unknown"
    assert first["n_classified"] == 1 and first["n_unknown"] == 1
    assert first["input_signature"]["mmseqs_legacy_rescue_sha256"] == hashlib.sha256(
        legacy_tsv.read_bytes()
    ).hexdigest()
    assert first["input_signature"]["mmseqs_legacy_rescue_algorithm_version"] == stage.EXTERNAL_RESCUE_ALGORITHM_VERSION
    assert "mmseqs_reference_rescue" not in first
    assert cache_state(case) == before
    assert case.calls() == []

    def unexpected(*args, **kwargs):
        pytest.fail("The refreshed matching legacy TSV checkpoint should be reused")

    monkeypatch.setattr(stage, "detect_complexity_from_fasta", unexpected)
    second = case.run()
    assert second["classified_library"] == first["classified_library"]
    assert second["n_classified"] == 1
    assert case.calls() == []


@pytest.mark.parametrize("change", ["changed", "removed", "emptied"])
def test_changed_or_removed_legacy_tsv_refreshes_classifications(many_query_case, change):
    case = many_query_case
    case.config.extra["classification"].pop("mmseqs_rescue_db")
    legacy_outputs(case)
    legacy_tsv = case.outdir / "mmseqs/mmseqs_rescue.tsv"
    legacy_tsv.write_text(
        "T-family-1#Unknown\tT-family-2#RC/Helitron\t0.95\t900\t0\t0\t1\t900\t1\t900\t1e-80\t1800\n"
    )
    first = case.run()
    assert read_table(first["classifications_tsv"])["T-family-1"]["order"] == "Helitron"
    before = cache_state(case)
    if change == "changed":
        legacy_tsv.write_text(legacy_tsv.read_text().replace("RC/Helitron", "DNA/hAT"))
    elif change == "removed":
        legacy_tsv.unlink()
    else:
        legacy_tsv.write_text("")

    refreshed = case.run()
    classified = read_table(refreshed["classifications_tsv"])
    assert classified["T-family-1"]["order"] == ("TIR" if change == "changed" else "Unknown")
    if change == "changed":
        assert refreshed["input_signature"]["mmseqs_legacy_rescue_sha256"] != first["input_signature"]["mmseqs_legacy_rescue_sha256"]
    else:
        assert "mmseqs_legacy_rescue_sha256" not in refreshed["input_signature"]
    assert cache_state(case) == before
    assert case.calls() == []


def test_legacy_tsv_changed_during_parse_cannot_publish_completion(many_query_case, monkeypatch):
    case = many_query_case
    case.config.extra["classification"].pop("mmseqs_rescue_db")
    legacy_outputs(case)
    legacy_tsv = case.outdir / "mmseqs/mmseqs_rescue.tsv"
    legacy_tsv.write_text(
        "T-family-1#Unknown\tT-family-2#RC/Helitron\t0.95\t900\t0\t0\t1\t900\t1\t900\t1e-80\t1800\n"
    )
    manifest = case.outdir / "classification.manifest.json"
    historical_manifest = manifest.read_bytes()
    parse = stage.parse_mmseqs_tsv

    def changed_after_parse(path, **kwargs):
        hits = parse(path, **kwargs)
        path.write_text(path.read_text().replace("RC/Helitron", "DNA/hAT"))
        return hits

    monkeypatch.setattr(stage, "parse_mmseqs_tsv", changed_after_parse)
    with pytest.raises(RuntimeError, match="Legacy MMseqs rescue TSV changed"):
        case.run()
    assert manifest.read_bytes() == historical_manifest
    assert case.calls() == []
