import csv
import json
import logging
from pathlib import Path

import pytest
import yaml

from drayte.pipeline import reclassify, reclassify_recovery as recovery
from drayte.pipeline.config import PipelineConfig
from drayte.pipeline.curation import write_final_curated_library
from drayte.storage.policy import output_lock


def make_run(tmp_path, families=3, reversals=1):
    root = tmp_path / "Mmex_DRayTEannotator2"
    sources = root / "extension/final_consensuses"
    sources.mkdir(parents=True)
    directory = root / "reclassify"
    directory.mkdir()
    merged = directory / "Mmex_extended_rep.fa"
    original_text = ""
    classified_text = ""
    rows = []
    for index in range(families):
        family = f"Mmex-fam{index:04d}"
        raw = f">{family}#Unknown\nAACG\n"
        (sources / f"{family}_rep.fa").write_text(raw)
        original_text += raw
        header = family + ("#DNA/Helitron" if index == 0 else "#LINE/CR1")
        classified_text += f">{header}\n" + ("CGTT" if index < reversals else "AACG") + "\n"
        rows.append({"name": family, "original_header": header, "consensus_len": 4})
    merged.write_text(original_text)
    classified = Path(str(merged) + ".classified")
    classified.write_text(classified_text)
    receipt = directory / "reclassify.manifest.json"
    receipt.write_text(json.dumps({"stage": "reclassify", "outdir": str(directory),
                                   "merged_library": str(merged), "classified_library": str(classified),
                                   "n_input_fastas": families}))
    curation = root / "curation"
    (curation / "prioritize").mkdir(parents=True)
    table = curation / "prioritize/family_table.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "original_header", "consensus_len"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    library = curation / "Final.RepeatModeler.Lib.fa"
    write_final_curated_library(classified, library, logging.getLogger("test.recovery"))
    curation_receipt = curation / "curation.manifest.json"
    curation_receipt.write_text(json.dumps({"stage": "curation", "outdir": str(curation),
                                            "final_library": str(library), "curated_library": str(library),
                                            "family_table": str(table)}))
    # Reproduce the user's already-damaged legacy chronology: the same merged
    # input is rewritten after both historical completion receipts.
    merged.write_text(original_text)
    classified.write_text(original_text)
    marker = directory / "repeatclassifier.checkpoint.json"
    marker.write_text(json.dumps({"schema_version": 1, "status": "running",
                                  "input_fasta": str(merged), "input_sha256": reclassify.sha256_file(merged),
                                  "output_fasta": str(classified), "repeatclassifier_bin": "RepeatClassifier"}))
    config = PipelineConfig(genome=str(tmp_path / "genome.fa"), outdir=str(root), species="Mmex")
    return config


def snapshot(root):
    return {str(path): (path.read_bytes(), path.stat().st_mtime_ns, path.stat().st_ctime_ns)
            for path in root.rglob("*") if path.is_file()}


def test_preview_is_read_only_and_verifies_damaged_legacy_chronology(tmp_path):
    config = make_run(tmp_path)
    before = snapshot(config.outdir_path)
    report = recovery.audit_recovery(config)
    assert report["eligible"], report
    assert report["sequence_comparison"] == {"identical_sequences": 2, "reverse_complement_sequences": 1}
    assert report["current_input"]["mtime_ns"] > report["historical_reclassify_receipt"]["mtime_ns"]
    assert snapshot(config.outdir_path) == before
    assert not (config.outdir_path / ".storage").exists()


def test_mmex_587_42_recovery_restores_headers_and_resumes_without_external_process(tmp_path, monkeypatch):
    config = make_run(tmp_path, families=587, reversals=42)
    root = config.outdir_path
    source = root / "curation/Final.RepeatModeler.Lib.fa"
    table = root / "curation/prioritize/family_table.tsv"
    source_before = snapshot(root / "curation")
    old_output = (root / "reclassify/Mmex_extended_rep.fa.classified").read_bytes()
    def unexpected(*args, **kwargs):
        raise AssertionError("Recovery and subsequent resume must not execute RepeatClassifier")
    monkeypatch.setattr(reclassify.subprocess, "Popen", unexpected)
    report = recovery.apply_recovery(config)
    assert report["applied"]
    assert not report["fresh_classifier_execution"]
    assert report["sequence_comparison"]["reverse_complement_sequences"] == 42
    assert snapshot(root / "curation") == source_before
    output = root / "reclassify/Mmex_extended_rep.fa.classified"
    assert output.read_text().startswith(">Mmex-fam0000#DNA/Helitron\nCGTT\n")
    assert source.read_text().startswith(">Mmex-fam0000#RC/Helitron\nCGTT\n")
    assert table.exists()
    backups = Path(report["backup_directory"])
    assert (backups / output.name).read_bytes() == old_output
    assert json.loads((backups / "repeatclassifier.checkpoint.json").read_text())["status"] == "running"
    assert (root / "reclassify/repeatclassifier.recovery.json").is_file()
    assert recovery.audit_recovery(config)["already_complete"]
    result = reclassify.run(config, {"extended_library_dir": str(root / "extension/final_consensuses")}, logging.getLogger("test"))
    assert Path(result["classified_library"]) == output


@pytest.mark.parametrize("change", ["source", "curated", "headers", "missing_receipt", "input_hash"])
def test_changed_or_incomplete_evidence_is_rejected_without_mutation(tmp_path, change):
    config = make_run(tmp_path)
    root = config.outdir_path
    if change == "source":
        (root / "extension/final_consensuses/Mmex-fam0000_rep.fa").write_text(">Mmex-fam0000#Unknown\nTTTT\n")
    elif change == "curated":
        path = root / "curation/Final.RepeatModeler.Lib.fa"
        path.write_text(path.read_text().replace("CGTT", "TTTT"))
    elif change == "headers":
        path = root / "curation/prioritize/family_table.tsv"
        path.write_text(path.read_text().replace("#DNA/Helitron", "#LTR/Gypsy"))
    elif change == "missing_receipt":
        (root / "curation/curation.manifest.json").unlink()
    else:
        path = root / "reclassify/repeatclassifier.checkpoint.json"
        record = json.loads(path.read_text())
        record["input_sha256"] = "wrong"
        path.write_text(json.dumps(record))
    before = snapshot(root)
    assert not recovery.audit_recovery(config)["eligible"]
    with pytest.raises((ValueError, OSError)):
        recovery.apply_recovery(config)
    assert {p: value for p, value in snapshot(root).items() if "/.storage/" not in p} == before


def test_modern_completed_checkpoint_cannot_be_overridden(tmp_path):
    config = make_run(tmp_path)
    root = config.outdir_path
    merged = root / "reclassify/Mmex_extended_rep.fa"
    output = Path(str(merged) + ".classified")
    path = root / "reclassify/repeatclassifier.checkpoint.json"
    path.write_text(json.dumps(reclassify.checkpoint_data(merged, output, "different-executable")))
    before = snapshot(root)
    assert not recovery.audit_recovery(config)["eligible"]
    with pytest.raises(ValueError, match="modern completed checkpoint"):
        recovery.apply_recovery(config)
    assert {p: value for p, value in snapshot(root).items() if "/.storage/" not in p} == before


@pytest.mark.parametrize("change", ["sequence", "classification"])
def test_semantic_mismatch_is_rejected_even_with_consistent_receipt_dates(tmp_path, change):
    config = make_run(tmp_path)
    root = config.outdir_path
    if change == "sequence":
        path = root / "curation/Final.RepeatModeler.Lib.fa"
        path.write_text(path.read_text().replace("CGTT", "TTTT"))
        expected = "differences other than orientation"
    else:
        path = root / "curation/prioritize/family_table.tsv"
        path.write_text(path.read_text().replace("#DNA/Helitron", "#LTR/Gypsy"))
        expected = "classification is inconsistent"
    # Isolate semantic validation from the timestamp guard: the otherwise
    # well-formed completion receipt is now newer than the bad source.
    receipt = root / "curation/curation.manifest.json"
    receipt.write_text(receipt.read_text())
    before = snapshot(root)
    report = recovery.audit_recovery(config)
    assert not report["eligible"]
    assert expected in report["errors"][0]
    assert snapshot(root) == before


def test_active_writer_lock_prevents_recovery(tmp_path):
    config = make_run(tmp_path)
    with output_lock(config.outdir_path):
        before = snapshot(config.outdir_path)
        report = recovery.audit_recovery(config)
        assert not report["eligible"]
        with pytest.raises(RuntimeError, match="holds the output lock"):
            recovery.apply_recovery(config)
        assert snapshot(config.outdir_path) == before


def test_apply_rolls_back_if_checkpoint_commit_fails(tmp_path, monkeypatch):
    config = make_run(tmp_path)
    root = config.outdir_path
    mutable = [root / "reclassify" / name for name in
               ("Mmex_extended_rep.fa.classified", "repeatclassifier.checkpoint.json", "reclassify.manifest.json")]
    before = {path: path.read_bytes() for path in mutable}
    real_write = recovery.write_json_atomic
    def fail_checkpoint(path, record):
        if path.name == "repeatclassifier.checkpoint.json":
            raise OSError("simulated quota failure")
        return real_write(path, record)
    monkeypatch.setattr(recovery, "write_json_atomic", fail_checkpoint)
    with pytest.raises(OSError, match="quota failure"):
        recovery.apply_recovery(config)
    assert all(path.read_bytes() == content for path, content in before.items())


def test_complete_checkpoint_is_published_only_after_recovery_provenance(tmp_path, monkeypatch):
    config = make_run(tmp_path)
    sidecar = config.outdir_path / "reclassify/repeatclassifier.recovery.json"
    real_write = recovery.write_json_atomic
    def checked_write(path, record):
        if path.name == "repeatclassifier.checkpoint.json" and record["status"] == "complete":
            assert sidecar.is_file()
            assert json.loads(sidecar.read_text())["fresh_classifier_execution"] is False
        return real_write(path, record)
    monkeypatch.setattr(recovery, "write_json_atomic", checked_write)
    assert recovery.apply_recovery(config)["applied"]


def test_cli_preview_returns_failure_for_missing_evidence(tmp_path, capsys):
    config = make_run(tmp_path)
    (config.outdir_path / "curation/Final.RepeatModeler.Lib.fa").unlink()
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({"genome": config.genome, "outdir": config.outdir, "species": config.species}))
    assert recovery.main(["--config", str(path)]) == 1
    assert not json.loads(capsys.readouterr().out)["eligible"]
