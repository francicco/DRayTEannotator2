import copy
import csv
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from drayte.classification import run
from drayte.classification.mmseqs import RescueEvidence
from drayte.classification.models import Family


def evidence(family_id="query"):
    return RescueEvidence(
        family_id=family_id, rescue_class="Class_II", rescue_order="Helitron",
        rescue_superfamily="Helitron", rescue_identity=0.95, rescue_aln_len=1100,
        rescue_query_cov=0.91, rescue_target_cov=0.91, rescue_bits=2000,
        rescue_target="reference-Helitron",
    )


def family(**kwargs):
    base = dict(family_id="query", consensus_len=1200, n_copies=0,
                rescue_class="Class_II", rescue_order="Helitron",
                rescue_superfamily="Helitron", rescue_identity=0.95,
                rescue_aln_len=1100, rescue_bits=2000,
                rescue_target="reference-Helitron")
    return Family(**{**base, **kwargs})


def test_rescue_impact_counts_actual_label_changes_without_mutating_inputs():
    supported = family()
    before = copy.deepcopy(supported)
    results = run.classify_families([supported])
    assert results[0]["order"] == "Helitron"
    assert run.summarize_rescue_impact([supported], results) == {
        "families_with_rescue_evidence": 1,
        "classifications_changed_by_rescue": 1,
        "previously_unknown_rescued": 1,
    }
    assert supported == before


def test_rescue_impact_does_not_claim_a_label_change_when_domains_override():
    supported = family(domains={"RT"})
    results = run.classify_families([supported])
    assert results[0]["order"] == "LINE"
    assert run.summarize_rescue_impact([supported], results) == {
        "families_with_rescue_evidence": 1,
        "classifications_changed_by_rescue": 0,
        "previously_unknown_rescued": 0,
    }


@pytest.mark.parametrize("mode", [None, "one-way", "reciprocal"])
def test_standalone_external_reference_uses_reference_api_despite_old_tsv(
    tmp_path, monkeypatch, capsys, mode,
):
    monkeypatch.chdir(tmp_path)
    query = tmp_path / "queries.fa"
    query.write_text(">query#Unknown\n" + "A" * 1200 + "\n")
    reference = tmp_path / "references.fa"
    reference.write_text(">reference-Helitron#RC/Helitron\n" + "A" * 1200 + "\n")
    outdir = tmp_path / "mmseqs"
    outdir.mkdir()
    (outdir / "mmseqs_rescue.tsv").write_text("legacy forward-only cache\n")
    calls = []

    def reference_rescue(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(
            forward_hits=[object()], reverse_hits=[object()], reciprocal_hits=[object()],
            evidence={"query": evidence()}, report={"cache_reused": False},
        )

    monkeypatch.setattr(run, "run_reference_mmseqs_rescue", reference_rescue)
    monkeypatch.setattr(run, "parse_mmseqs_tsv", lambda *args, **kwargs: pytest.fail("Unverified TSV reused"))
    output, detailed = tmp_path / "classified.tsv", tmp_path / "evidence.tsv"
    argv = [
        "drayte-classify", "--fasta", str(query), "--mmseqs-rescue-db", str(reference),
        "--mmseqs-outdir", str(outdir), "--output", str(output),
        "--evidence-output", str(detailed),
    ]
    if mode is not None:
        argv.extend(["--mmseqs-rescue-mode", mode])
    monkeypatch.setattr(sys, "argv", argv)
    run.main()
    assert len(calls) == 1
    assert Path(calls[0]["reference_fasta"]) == reference
    assert calls[0]["mode"] == (mode or "one-way")
    assert calls[0]["min_score_fraction"] == 0.95
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1 and rows[0]["family_id"] == "query"
    assert rows[0]["order"] == "Helitron"
    with detailed.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["rescue_target"] == "reference-Helitron"
    assert "families_with_rescue_evidence=1" in capsys.readouterr().err


@pytest.mark.parametrize("args", [
    ["--fasta", "query.fa", "--mmseqs-rescue-tsv", "supplied.tsv"],
    ["--input", "families.tsv"],
])
def test_standalone_rejects_ambiguous_or_sequence_free_reference_request(monkeypatch, args):
    monkeypatch.setattr(run, "run_reference_mmseqs_rescue", lambda **kwargs: pytest.fail("Unexpected search"))
    monkeypatch.setattr(sys, "argv", ["drayte-classify", *args,
                                     "--mmseqs-rescue-db", "reference.fa", "--output", "out.tsv"])
    with pytest.raises(SystemExit) as error:
        run.main()
    assert error.value.code == 2


def test_standalone_without_reference_keeps_legacy_tsv_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    query = tmp_path / "queries.fa"
    query.write_text(">query#Unknown\n" + "A" * 1200 + "\n")
    cached = tmp_path / "legacy.tsv"
    cached.write_text("legacy cache\n")
    calls = []
    monkeypatch.setattr(run, "run_reference_mmseqs_rescue", lambda **kwargs: pytest.fail("Unexpected reference search"))
    monkeypatch.setattr(run, "parse_mmseqs_tsv", lambda path, **kwargs: calls.append(path) or [])
    monkeypatch.setattr(sys, "argv", ["drayte-classify", "--fasta", str(query),
                                     "--mmseqs-rescue-tsv", str(cached),
                                     "--output", str(tmp_path / "out.tsv")])
    run.main()
    assert calls == [cached]


@pytest.mark.parametrize("fraction", ["0", "-0.1", "1.1", "nan", "inf"])
def test_standalone_rejects_invalid_reference_score_fraction(monkeypatch, fraction):
    monkeypatch.setattr(run, "run_reference_mmseqs_rescue", lambda **kwargs: pytest.fail("Unexpected search"))
    monkeypatch.setattr(sys, "argv", [
        "drayte-classify", "--fasta", "query.fa", "--mmseqs-rescue-db", "reference.fa",
        "--mmseqs-rescue-score-fraction", fraction, "--output", "out.tsv",
    ])
    with pytest.raises(SystemExit) as error:
        run.main()
    assert error.value.code == 2


def test_standalone_passes_configured_score_fraction_and_rescues_multiple_queries(
    tmp_path, monkeypatch,
):
    monkeypatch.chdir(tmp_path)
    query = tmp_path / "queries.fa"
    query.write_text(">query#Unknown\n" + "A" * 1200 + "\n>query-2#Unknown\n" + "C" * 1200 + "\n")
    reference = tmp_path / "references.fa"
    reference.write_text(">reference-Helitron#RC/Helitron\n" + "A" * 1200 + "\n")
    calls = []

    def reference_rescue(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(
            forward_hits=[object(), object()], reverse_hits=[], reciprocal_hits=[],
            evidence={name: evidence(name) for name in ("query", "query-2")},
            report={"cache_reused": False},
        )

    monkeypatch.setattr(run, "run_reference_mmseqs_rescue", reference_rescue)
    output = tmp_path / "classifications.tsv"
    monkeypatch.setattr(sys, "argv", [
        "drayte-classify", "--fasta", str(query), "--mmseqs-rescue-db", str(reference),
        "--mmseqs-rescue-score-fraction", "0.9", "--output", str(output),
    ])
    run.main()
    assert calls[0]["mode"] == "one-way"
    assert calls[0]["min_score_fraction"] == 0.9
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["family_id"] for row in rows} == {"query", "query-2"}
    assert all(row["order"] == "Helitron" for row in rows)
