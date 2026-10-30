import json
import logging
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from drayte.pipeline import reclassify


def legacy_files(tmp_path, *, transformed=True):
    merged = tmp_path / "merged.fa"
    output = tmp_path / "merged.fa.classified"
    receipt = tmp_path / "reclassify.manifest.json"
    merged.write_text(">fam1#Unknown\nACGT\n")
    output.write_text(">fam1#DNA\n" + ("AAAA" if transformed else "ACGT") + "\n")
    receipt.write_text(json.dumps({"stage": "reclassify", "merged_library": str(merged),
                                   "classified_library": str(output)}))
    for index, path in enumerate((merged, output, receipt), 1):
        os.utime(path, ns=(index * 1_000_000_000, index * 1_000_000_000))
    return merged, output, receipt


def fake_success(monkeypatch, output):
    calls = []
    class Process:
        stdout = []
        def wait(self):
            output.write_text(">fam1#DNA\nACGT\n")
            return 0
    def popen(*args, **kwargs):
        calls.append(args)
        return Process()
    monkeypatch.setattr(reclassify.subprocess, "Popen", popen)
    return calls


def test_transformed_legacy_output_with_success_receipt_is_migrated(tmp_path, monkeypatch):
    merged, output, _ = legacy_files(tmp_path)
    def unexpected(*args, **kwargs):
        raise AssertionError("Valid transformed legacy output must not rerun")
    monkeypatch.setattr(reclassify.subprocess, "Popen", unexpected)
    original = output.read_bytes()
    assert reclassify.run_repeatclassifier(merged, tmp_path, "RepeatClassifier", logging.getLogger("test")) == output
    assert output.read_bytes() == original
    assert not list(tmp_path.glob("*.stale*"))
    assert reclassify.checkpoint_is_valid(tmp_path / "repeatclassifier.checkpoint.json", merged, output, "RepeatClassifier")


def test_changed_executable_cannot_fall_back_to_legacy_receipt(tmp_path, monkeypatch):
    merged, output, _ = legacy_files(tmp_path, transformed=False)
    (tmp_path / "repeatclassifier.checkpoint.json").write_text(json.dumps(
        reclassify.checkpoint_data(merged, output, "old-RepeatClassifier")))
    calls = fake_success(monkeypatch, output)
    reclassify.run_repeatclassifier(merged, tmp_path, "new-RepeatClassifier", logging.getLogger("test"))
    assert len(calls) == 1


@pytest.mark.parametrize("condition", ["no_receipt", "newer_input", "newer_output", "running_marker", "invalid_marker"])
def test_unproven_or_unfinished_outputs_are_not_migrated(tmp_path, monkeypatch, condition):
    merged, output, receipt = legacy_files(tmp_path)
    if condition == "no_receipt":
        receipt.unlink()
    elif condition == "newer_input":
        os.utime(merged, ns=(4_000_000_000, 4_000_000_000))
    elif condition == "newer_output":
        os.utime(output, ns=(4_000_000_000, 4_000_000_000))
    elif condition == "running_marker":
        (tmp_path / "repeatclassifier.checkpoint.json").write_text('{"status": "running"}')
    else:
        (tmp_path / "repeatclassifier.checkpoint.json").write_text("broken JSON")
    calls = fake_success(monkeypatch, output)
    reclassify.run_repeatclassifier(merged, tmp_path, "RepeatClassifier", logging.getLogger("test"))
    assert len(calls) == 1


def test_unchanged_merge_preserves_bytes_and_timestamp(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "fam1.fa").write_text(">fam1#Unknown\nACGT\n")
    merged = tmp_path / "merged.fa"
    reclassify.gather_extended_consensi(source, "fam", merged, logging.getLogger("test"))
    os.utime(merged, ns=(1234567890123456789, 1234567890123456789))
    original = merged.read_bytes(), merged.stat().st_mtime_ns
    reclassify.gather_extended_consensi(source, "fam", merged, logging.getLogger("test"))
    assert (merged.read_bytes(), merged.stat().st_mtime_ns) == original


def test_mmex_587_families_42_transformations_resume_without_rerun(tmp_path, monkeypatch):
    root = tmp_path / "run"
    source = root / "extension/final_consensuses"
    source.mkdir(parents=True)
    outdir = root / "reclassify"
    outdir.mkdir()
    for index in range(587):
        (source / f"Mmex-fam{index:04d}_rep.fa").write_text(f">Mmex-fam{index:04d}#Unknown\nACGT\n")
    merged = outdir / "Mmex_extended_rep.fa"
    reclassify.gather_extended_consensi(source, "Mmex", merged, logging.getLogger("test"))
    output = outdir / "Mmex_extended_rep.fa.classified"
    output.write_text("".join(f">Mmex-fam{index:04d}#DNA\n" + ("AAAA" if index < 42 else "ACGT") + "\n"
                              for index in range(587)))
    receipt = outdir / "reclassify.manifest.json"
    receipt.write_text(json.dumps({"stage": "reclassify", "merged_library": str(merged),
                                   "classified_library": str(output)}))
    for index, path in enumerate((merged, output, receipt), 1):
        os.utime(path, ns=(index * 1_000_000_000, index * 1_000_000_000))
    def unexpected(*args, **kwargs):
        raise AssertionError("The completed 587-family legacy result should resume")
    monkeypatch.setattr(reclassify.subprocess, "Popen", unexpected)
    result = reclassify.run(SimpleNamespace(outdir_path=root, species="Mmex", extra={}),
                            {"extended_library_dir": str(source)}, logging.getLogger("test"))
    assert Path(result["classified_library"]) == output
    assert reclassify.compare_classified_output(merged, output)["changed_sequences"] == 42


@pytest.mark.parametrize("content", ["", ">fam2#DNA\nACGT\n", ">fam1#DNA\nACG\n"])
def test_invalid_legacy_family_set_cannot_be_adopted(tmp_path, monkeypatch, content):
    merged, output, _ = legacy_files(tmp_path)
    output.write_text(content)
    os.utime(output, ns=(2_000_000_000, 2_000_000_000))
    calls = fake_success(monkeypatch, output)
    reclassify.run_repeatclassifier(merged, tmp_path, "RepeatClassifier", logging.getLogger("test"))
    assert len(calls) == 1


def test_legacy_hash_receipt_rejects_changed_input_even_if_timestamp_preserved(tmp_path, monkeypatch):
    merged, output, receipt = legacy_files(tmp_path)
    data = json.loads(receipt.read_text())
    data.update(input_sha256=reclassify.sha256_file(merged), output_sha256=reclassify.sha256_file(output))
    receipt.write_text(json.dumps(data))
    merged.write_text(">fam1#Unknown\nTTTT\n")
    os.utime(merged, ns=(1_000_000_000, 1_000_000_000))
    calls = fake_success(monkeypatch, output)
    reclassify.run_repeatclassifier(merged, tmp_path, "RepeatClassifier", logging.getLogger("test"))
    assert len(calls) == 1


def test_failed_attempt_keeps_running_marker_and_cannot_use_legacy_receipt(tmp_path, monkeypatch):
    merged, output, _ = legacy_files(tmp_path, transformed=False)
    checkpoint = tmp_path / "repeatclassifier.checkpoint.json"
    checkpoint.write_text('{"status": "running"}')
    calls = []
    class Process:
        stdout = []
        def wait(self):
            output.write_text(">fam1#DNA\nACGT\n")
            # Even if output metadata resembles a completed legacy run, the
            # explicit running marker must prevent fallback after failure.
            os.utime(output, ns=(2_000_000_000, 2_000_000_000))
            return 1 if len(calls) == 1 else 0
    def popen(*args, **kwargs):
        calls.append(args)
        return Process()
    monkeypatch.setattr(reclassify.subprocess, "Popen", popen)
    with pytest.raises(RuntimeError, match="exit code 1"):
        reclassify.run_repeatclassifier(merged, tmp_path, "RepeatClassifier", logging.getLogger("test"))
    assert json.loads(checkpoint.read_text())["status"] == "running"
    reclassify.run_repeatclassifier(merged, tmp_path, "RepeatClassifier", logging.getLogger("test"))
    assert len(calls) == 2
    assert json.loads(checkpoint.read_text())["status"] == "complete"


def test_interrupted_merge_keeps_previous_library(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    fasta = source / "fam1.fa"
    fasta.write_text(">fam1\nACGT\n")
    merged = tmp_path / "merged.fa"
    reclassify.gather_extended_consensi(source, "fam", merged, logging.getLogger("test"))
    original = merged.read_bytes(), merged.stat().st_mtime_ns
    fasta.write_text(">fam1\nTTTT\n")
    def fail(*args, **kwargs):
        raise OSError("simulated write failure")
    monkeypatch.setattr(reclassify.shutil, "copyfileobj", fail)
    with pytest.raises(OSError, match="write failure"):
        reclassify.gather_extended_consensi(source, "fam", merged, logging.getLogger("test"))
    assert (merged.read_bytes(), merged.stat().st_mtime_ns) == original
    assert not list(tmp_path.glob(".merged.fa.*.partial"))
