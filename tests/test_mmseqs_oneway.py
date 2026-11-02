"""Curated reference classification permits many query models per reference."""
import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from drayte.classification import mmseqs


def row(query, target, *, bits=1800, identity=.91, length=900):
    return f"{query}\t{target}\t{identity}\t{length}\t0\t0\t1\t{length}\t1\t{length}\t1e-50\t{bits}\n"


@pytest.fixture
def oneway(tmp_path):
    query = tmp_path / "query.fa"
    reference = tmp_path / "reference.fa"
    query.write_text(">q1#Unknown\n" + "ACGT" * 250 + "\n")
    reference.write_text(">r1#RC/Helitron\n" + "ACGT" * 250 + "\n")
    config = tmp_path / "fake.json"
    calls = tmp_path / "calls.jsonl"
    binary = tmp_path / "fake-mmseqs"
    binary.write_text(
        f"#!{sys.executable}\n"
        "import json, sys\nfrom pathlib import Path\n"
        f"config = json.loads(Path({str(config)!r}).read_text())\n"
        f"with Path({str(calls)!r}).open('a') as log: log.write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "source = Path(sys.argv[2]).name\n"
        "if config.get('fail') == source: sys.exit(7)\n"
        "Path(sys.argv[4]).write_text(config.get(source, ''))\n"
    )
    binary.chmod(0o755)
    directory = tmp_path / "rescue"

    def settings(forward=None, fail=None):
        config.write_text(json.dumps({
            "query.fa": row("q1#Unknown", "r1#RC/Helitron") if forward is None else forward,
            "fail": fail,
        }))

    def run(**kwargs):
        return mmseqs.run_reference_mmseqs_rescue(query, reference, directory, str(binary), **kwargs)

    settings()
    return SimpleNamespace(query=query, reference=reference, config=config,
                           calls=calls, directory=directory, run=run, settings=settings)


def add_reference(case, name):
    with case.reference.open("a") as handle:
        handle.write(f">{name}\n" + "ACGT" * 250 + "\n")


def test_many_query_models_receive_same_reference_label_in_one_search(oneway):
    case = oneway
    with case.query.open("a") as handle:
        handle.write(">q2#Unknown\n" + "ACGT" * 250 + "\n")
    case.settings(forward=row("q1", "r1") + row("q2", "r1", bits=1700), fail="reference.fa")
    result = case.run()
    assert set(result.evidence) == {"q1", "q2"}
    assert all(item.rescue_order == "Helitron" for item in result.evidence.values())
    assert len(case.calls.read_text().splitlines()) == 1
    assert result.report["mode"] == "one-way"
    assert not any("reverse" in key or "reciprocal" in key for key in result.report)
    assert not (case.directory / "mmseqs_rescue.reverse.tsv").exists()
    assert result.report["query_decisions"]["q2"]["status"] == "rescued"


def test_identical_clean_ids_across_libraries_are_valid(oneway):
    oneway.reference.write_text(">q1#RC/Helitron\n" + "ACGT" * 250 + "\n")
    oneway.settings(forward=row("q1#Unknown", "q1#RC/Helitron"))
    assert oneway.run().evidence["q1"].rescue_order == "Helitron"


def test_same_label_ties_are_accepted_deterministically(oneway):
    add_reference(oneway, "r2#RC/Helitron")
    oneway.settings(forward=row("q1", "r2") + row("q1", "r1"))
    result = oneway.run()
    assert result.evidence["q1"].rescue_target == "r1#RC/Helitron"
    assert result.evidence["q1"].rescue_superfamily == "Helitron"
    assert len(result.report["query_decisions"]["q1"]["near_best_targets"]) == 2


def test_conflicting_near_best_orders_abstain(oneway):
    add_reference(oneway, "r2#DNA/hAT")
    oneway.settings(forward=row("q1", "r1", bits=2000) + row("q1", "r2", bits=1950))
    result = oneway.run()
    assert not result.evidence
    assert result.report["query_decisions"]["q1"]["reason"] == "conflicting_class_or_order"


def test_same_order_conflicting_superfamilies_downgrades_assignment(oneway):
    oneway.reference.write_text(">r1#DNA/hAT\n" + "ACGT" * 250 + "\n")
    add_reference(oneway, "r2#DNA/TcMar-Mariner")
    oneway.settings(forward=row("q1", "r1") + row("q1", "r2", bits=1750))
    result = oneway.run()
    assert result.evidence["q1"].rescue_order == "TIR"
    assert result.evidence["q1"].rescue_superfamily == "Unknown"


def test_weaker_conflicting_target_outside_score_window_is_ignored(oneway):
    add_reference(oneway, "r2#DNA/hAT")
    oneway.settings(forward=row("q1", "r1", bits=2000) + row("q1", "r2", bits=1899))
    assert oneway.run().evidence["q1"].rescue_order == "Helitron"


@pytest.mark.parametrize("label", ["Unknown", "Unrecognized", ""])
def test_unclassified_competing_reference_cannot_be_skipped(oneway, label):
    add_reference(oneway, f"r2#{label}" if label else "r2")
    oneway.settings(forward=row("q1", "r2", bits=2000) + row("q1", "r1", bits=1950))
    result = oneway.run()
    assert not result.evidence
    assert result.report["query_decisions"]["q1"]["reason"] == "unclassified_competing_reference"


def test_duplicate_hsps_are_not_summed_to_outvote_stronger_reference(oneway):
    add_reference(oneway, "r2#DNA/hAT")
    oneway.settings(forward=row("q1", "r1", bits=2000)
                    + row("q1", "r2", bits=1200) + row("q1", "r2", bits=1200))
    result = oneway.run()
    assert result.evidence["q1"].rescue_order == "Helitron"
    assert result.report["query_decisions"]["q1"]["qualified_targets"] == 2


@pytest.mark.parametrize("options,failed", [
    ({"identity": .79}, "identity"),
    ({"length": 499}, "alignment_length"),
    ({"length": 600}, "query_coverage"),
])
def test_retained_filters_explain_rejection(oneway, options, failed):
    oneway.settings(forward=row("q1", "r1", **options))
    result = oneway.run(min_query_cov=.8)
    assert not result.evidence
    decision = result.report["query_decisions"]["q1"]
    assert decision["reason"] == "no_qualified_hit"
    assert failed in decision["rejected_hits"][0]["failed_filters"]


def test_shorter_target_coverage_keeps_order_without_superfamily(oneway):
    oneway.reference.write_text(">r1#RC/Helitron\n" + "ACGT" * 500 + "\n")
    result = oneway.run()
    assert result.evidence["q1"].rescue_order == "Helitron"
    assert result.evidence["q1"].rescue_superfamily == "Unknown"


def test_reference_coverage_failure_is_reported(oneway):
    oneway.reference.write_text(">r1#RC/Helitron\n" + "ACGT" * 1000 + "\n")
    result = oneway.run()
    assert not result.evidence
    assert result.report["query_decisions"]["q1"]["rejected_hits"][0]["failed_filters"] == ["reference_coverage"]


@pytest.mark.parametrize("bits", [0, -10])
def test_nonpositive_bit_scores_cannot_support_rescue(oneway, bits):
    oneway.settings(forward=row("q1", "r1", bits=bits))
    result = oneway.run()
    assert not result.evidence
    decision = result.report["query_decisions"]["q1"]
    assert decision["reason"] == "no_positive_bit_score"
    assert decision["rejected_hits"][0]["failed_filters"] == ["nonpositive_bit_score"]


def test_score_window_change_is_signed_and_changes_ambiguity(oneway):
    add_reference(oneway, "r2#DNA/hAT")
    oneway.settings(forward=row("q1", "r1", bits=2000) + row("q1", "r2", bits=1950))
    assert not oneway.run().evidence
    result = oneway.run(min_score_fraction=.99)
    assert result.evidence["q1"].rescue_order == "Helitron"
    assert not result.report["cache_reused"]
    assert len(oneway.calls.read_text().splitlines()) == 2


def test_historical_reciprocal_cache_is_rebuilt_and_old_reverse_is_untouched(oneway):
    first = oneway.run()
    checkpoint_path = Path(first.report["checkpoint"])
    checkpoint = json.loads(checkpoint_path.read_text())
    checkpoint["input_signature"]["mode"] = "reciprocal"
    checkpoint["reverse_sha256"] = "historical digest"
    checkpoint_path.write_text(json.dumps(checkpoint))
    old_reverse = oneway.directory / "mmseqs_rescue.reverse.tsv"
    old_reverse.write_text("historical reverse content\n")
    before = old_reverse.read_bytes(), old_reverse.stat().st_mtime_ns
    rebuilt = oneway.run()
    assert not rebuilt.report["cache_reused"]
    assert len(oneway.calls.read_text().splitlines()) == 2
    assert (old_reverse.read_bytes(), old_reverse.stat().st_mtime_ns) == before
    assert not any("reverse" in key for key in rebuilt.report)
    cached = oneway.run(threads=8)
    assert cached.report["cache_reused"]
    assert len(oneway.calls.read_text().splitlines()) == 2


def test_completed_v2_oneway_cache_with_old_null_reverse_field_is_reused(oneway):
    first = oneway.run()
    checkpoint_path = Path(first.report["checkpoint"])
    checkpoint = json.loads(checkpoint_path.read_text())
    checkpoint["reverse_sha256"] = None
    checkpoint_path.write_text(json.dumps(checkpoint))
    # This legacy path may even be a dangling symlink. It is never opened or
    # validated by the one-way cache.
    old_reverse = oneway.directory / "mmseqs_rescue.reverse.tsv"
    old_reverse.symlink_to(oneway.directory / "missing-historical-reverse")
    cached = oneway.run()
    assert cached.report["cache_reused"]
    assert len(oneway.calls.read_text().splitlines()) == 1
    assert old_reverse.is_symlink()
    assert not any("reverse" in key for key in cached.report)


def test_empty_oneway_search_is_reusable_and_explains_missing_hits(oneway):
    oneway.settings(forward="")
    first = oneway.run()
    second = oneway.run(threads=8)
    assert not first.evidence and not second.evidence
    assert second.report["cache_reused"]
    assert first.report["query_decisions"]["q1"]["reason"] == "no_search_hit"
    assert len(oneway.calls.read_text().splitlines()) == 1


def test_oneway_publication_failure_restores_saved_outputs(oneway, monkeypatch):
    first = oneway.run()
    paths = [Path(first.report[key]) for key in ("forward_tsv", "report_path")]
    before = {path: path.read_bytes() for path in paths}
    replace = os.replace

    def fail_report(source, destination):
        if Path(destination) == paths[1] and Path(source).suffix == ".partial":
            raise OSError("simulated report publication failure")
        return replace(source, destination)

    monkeypatch.setattr(mmseqs.os, "replace", fail_report)
    with pytest.raises(OSError, match="publication failure"):
        oneway.run(min_identity=.81)
    assert all(path.read_bytes() == data for path, data in before.items())
    assert json.loads(Path(first.report["checkpoint"]).read_text())["status"] == "failed"


@pytest.mark.parametrize("fraction", [0, -1, 1.01, float("nan"), float("inf")])
def test_invalid_score_fraction_fails_before_search(oneway, fraction):
    with pytest.raises(ValueError, match="min_score_fraction"):
        oneway.run(min_score_fraction=fraction)
    assert not oneway.calls.exists()
