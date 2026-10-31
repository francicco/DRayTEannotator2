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


def test_standalone_external_reference_uses_reciprocal_api_despite_old_tsv(
    tmp_path, monkeypatch, capsys,
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

    def reciprocal(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(
            forward_hits=[object()], reverse_hits=[object()], reciprocal_hits=[object()],
            evidence={"query": evidence()}, report={"cache_reused": False},
        )

    monkeypatch.setattr(run, "run_reciprocal_mmseqs_rescue", reciprocal)
    monkeypatch.setattr(run, "parse_mmseqs_tsv", lambda *args, **kwargs: pytest.fail("Unverified TSV reused"))
    output, detailed = tmp_path / "classified.tsv", tmp_path / "evidence.tsv"
    monkeypatch.setattr(sys, "argv", [
        "drayte-classify", "--fasta", str(query), "--mmseqs-rescue-db", str(reference),
        "--mmseqs-outdir", str(outdir), "--output", str(output),
        "--evidence-output", str(detailed),
    ])
    run.main()
    assert len(calls) == 1
    assert Path(calls[0]["reference_fasta"]) == reference
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
    monkeypatch.setattr(run, "run_reciprocal_mmseqs_rescue", lambda **kwargs: pytest.fail("Unexpected search"))
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
    monkeypatch.setattr(run, "run_reciprocal_mmseqs_rescue", lambda **kwargs: pytest.fail("Unexpected reference search"))
    monkeypatch.setattr(run, "parse_mmseqs_tsv", lambda path, **kwargs: calls.append(path) or [])
    monkeypatch.setattr(sys, "argv", ["drayte-classify", "--fasta", str(query),
                                     "--mmseqs-rescue-tsv", str(cached),
                                     "--output", str(tmp_path / "out.tsv")])
    run.main()
    assert calls == [cached]
