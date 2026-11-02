"""Signed one-way reference searches, namespaces, and safe cache publication."""
import fcntl
import json
import os
import sys
from pathlib import Path

import pytest

from drayte.classification import mmseqs


def row(query, target, bits=1800, identity=.91, length=900):
    return f"{query}\t{target}\t{identity}\t{length}\t0\t0\t1\t{length}\t1\t{length}\t1e-50\t{bits}\n"


@pytest.fixture
def search(tmp_path):
    query = tmp_path / "query.fa"
    reference = tmp_path / "reference.fa"
    query.write_text(">query-family#Unknown\n" + "ACGT" * 250 + "\n")
    reference.write_text(">reference-family#RC/Helitron\n" + "ACGT" * 250 + "\n")
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

    def settings(forward=None, fail=None):
        config.write_text(json.dumps({
            "query.fa": row("query-family#Unknown", "reference-family#RC/Helitron") if forward is None else forward,
            "fail": fail,
        }))

    def run(**kwargs):
        return mmseqs.run_reference_mmseqs_rescue(query, reference, tmp_path / "rescue", str(binary), **kwargs)

    settings()
    return type("Search", (), {"query": query, "reference": reference, "config": config,
                               "calls": calls, "binary": binary, "directory": tmp_path / "rescue",
                               "settings": staticmethod(settings), "run": staticmethod(run)})()


def test_one_search_rescues_disjoint_ids_without_reference_families(search):
    result = search.run()
    calls = [json.loads(line) for line in search.calls.read_text().splitlines()]
    assert len(calls) == 1
    assert calls[0][1:3] == [str(search.query), str(search.reference)]
    assert calls[0][calls[0].index("--format-output") + 1] == mmseqs.MMSEQS_OUTPUT_FORMAT
    assert set(result.evidence) == {"query-family"}
    assert result.evidence["query-family"].rescue_order == "Helitron"
    assert result.evidence["query-family"].rescue_superfamily == "Helitron"
    assert result.report["families_with_rescue"] == result.report["forward_hits"] == 1
    assert not any("reverse" in key or "reciprocal" in key for key in result.report)
    checkpoint = json.loads(Path(result.report["checkpoint"]).read_text())
    assert checkpoint["status"] == "complete"
    assert not any("reverse" in key for key in checkpoint)
    assert "query_to_reference" in (search.directory / "mmseqs_rescue.log").read_text()
    assert "reference_to_query" not in (search.directory / "mmseqs_rescue.log").read_text()


def test_namespaces_allow_same_clean_id_in_different_libraries(search):
    search.reference.write_text(">query-family#RC/Helitron\n" + "ACGT" * 250 + "\n")
    search.settings(forward=row("query-family#Unknown", "query-family#RC/Helitron"))
    result = search.run()
    assert set(result.evidence) == {"query-family"}
    assert result.evidence["query-family"].rescue_order == "Helitron"


def test_classification_header_is_taken_from_signed_reference_not_tool_output(search):
    search.settings(forward=row("query-family", "reference-family"))
    result = search.run()
    assert result.evidence["query-family"].rescue_target == "reference-family#RC/Helitron"


def test_query_and_reference_coverage_use_their_own_lengths(search):
    search.reference.write_text(">reference-family#RC/Helitron\n" + "ACGT" * 375 + "\n")
    result = search.run(min_query_cov=.8, min_target_cov=.5)
    assert result.forward_hits[0].query_cov == .9
    assert result.forward_hits[0].target_cov == .6
    assert result.evidence
    assert not search.run(min_target_cov=.7).evidence


def test_forward_best_support_and_label_conflicts_control_rescue(search):
    search.reference.write_text(search.reference.read_text() + ">second-reference#DNA/hAT\n" + "ACGT" * 250 + "\n")
    search.settings(forward=row("query-family#Unknown", "reference-family#RC/Helitron", bits=1000)
                    + row("query-family#Unknown", "second-reference#DNA/hAT", bits=2000))
    assert search.run().evidence["query-family"].rescue_order == "TIR"
    search.settings(forward=row("query-family#Unknown", "reference-family#RC/Helitron", bits=2000)
                    + row("query-family#Unknown", "second-reference#DNA/hAT", bits=2000))
    assert not search.run(min_identity=.81).evidence


def test_zero_hit_outputs_are_complete_reusable_cache(search):
    search.settings(forward="")
    first = search.run()
    second = search.run(threads=12)
    assert not first.evidence and not second.evidence
    assert first.report["cache_reused"] is False
    assert second.report["cache_reused"] is True
    assert len(search.calls.read_text().splitlines()) == 1


@pytest.mark.parametrize("changed", ["query", "reference", "filter", "corrupt_output"])
def test_input_or_filter_or_output_change_invalidates_cache(search, changed):
    first = search.run()
    options = {}
    if changed in {"query", "reference"}:
        path = getattr(search, changed)
        path.write_text(path.read_text().replace("ACGT", "ACGA", 1))
    elif changed == "filter":
        options = {"min_identity": .95}
    else:
        Path(first.report["forward_tsv"]).write_text("untrusted output\n")
    second = search.run(**options)
    assert not second.report["cache_reused"]
    assert len(search.calls.read_text().splitlines()) == 2


def test_unsigned_one_direction_file_is_never_adopted(search):
    search.directory.mkdir()
    (search.directory / "mmseqs_rescue.tsv").write_text(row("other", "unrelated"))
    assert search.run().evidence
    assert len(search.calls.read_text().splitlines()) == 1


def test_failed_search_keeps_previous_outputs_and_incomplete_marker(search):
    first = search.run()
    files = [Path(first.report[key]) for key in ("forward_tsv", "report_path")]
    before = {path: path.read_bytes() for path in files}
    search.reference.write_text(search.reference.read_text().replace("ACGT", "ACGA", 1))
    search.settings(fail="query.fa")
    with pytest.raises(Exception):
        search.run()
    assert all(path.read_bytes() == data for path, data in before.items())
    checkpoint = json.loads(Path(first.report["checkpoint"]).read_text())
    assert checkpoint["status"] == "failed"
    search.settings()
    assert search.run().evidence
    assert len(search.calls.read_text().splitlines()) == 3


def test_publication_failure_rolls_back_both_existing_outputs(search, monkeypatch):
    first = search.run()
    files = [Path(first.report[key]) for key in ("forward_tsv", "report_path")]
    before = {path: path.read_bytes() for path in files}
    replace = os.replace

    def fail_report(source, destination):
        if Path(destination) == files[1] and Path(source).suffix == ".partial":
            raise OSError("simulated publication failure")
        return replace(source, destination)

    monkeypatch.setattr(mmseqs.os, "replace", fail_report)
    with pytest.raises(OSError, match="publication failure"):
        search.run(min_identity=.81)
    assert all(path.read_bytes() == data for path, data in before.items())
    assert json.loads(Path(first.report["checkpoint"]).read_text())["status"] == "failed"


def test_duplicate_clean_ids_fail_before_search(search):
    search.reference.write_text(search.reference.read_text() + ">reference-family#DNA/hAT\nAAAA\n")
    with pytest.raises(ValueError, match="Duplicate"):
        search.run()
    assert not search.calls.exists()


def test_unknown_output_ids_cannot_become_query_evidence(search):
    search.settings(forward=row("reference-family#RC/Helitron", "reference-family#RC/Helitron"))
    with pytest.raises(ValueError, match="unexpected sequence ID"):
        search.run()
    assert json.loads((search.directory / "mmseqs_rescue.checkpoint.json").read_text())["status"] == "failed"


def test_concurrent_rescue_lock_blocks_before_tools_run(search):
    search.directory.mkdir()
    with (search.directory / "mmseqs_rescue.lock").open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(RuntimeError, match="holds the output lock"):
            search.run()
    assert not search.calls.exists()
