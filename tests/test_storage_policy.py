import json
import logging
import sys
from pathlib import Path

import pytest

from drayte.pipeline import classification, reclassify, runner
from drayte.pipeline.config import PipelineConfig
from drayte.storage.archive import ArchiveStore
from drayte.storage.policy import (
    assert_safe_deletions, cleanup_stage_scratch, compact_completed_extension, compact_run, output_lock,
    plan_compaction, storage_mode, validate_completed_run,
)
from drayte.utils.repeatmasker_checkpoint import input_signature, write_checkpoint


@pytest.fixture
def completed_run(tmp_path):
    root = tmp_path / "Test_DRayTEannotator2"
    config = PipelineConfig(genome=str(tmp_path / "genome.fa"), outdir=str(root), species="Test",
                            extra={"storage": {"mode": "compact"}, "family_inspection": {"enabled": False}})

    def write(relative, text):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    sequence = ">fam1#Unknown\nACGT\n"
    genome = write("discovery/assemblies_dir/Test.fa", ">chr1\nACGTACGT\n")
    Path(config.genome).write_text(genome.read_text())
    raw = write("discovery/rmodeler_dir/consensi.fa.classified", sequence)
    write("discovery/discovery.manifest.json", json.dumps({"raw_library": str(raw)}))
    write("discovery/rmasker_dir/Test.fa.out", "repeat hits\n")
    write("discovery/rmasker_dir/Test.fa.tbl", "repeat table\n")
    write("discovery/rmasker_dir/Test.fa.masked", "ACGT" * 2000)
    write("discovery/rmodeler_dir/RM_old/round-1/work.fa", "unused" * 1000)
    write("extension/blastfiles/Test-families.mod.fa", sequence)
    write("extension/extensionwork/fam1/repam-ranges.tsv", "1\t4\n")
    write("extension/extensionwork/fam1/repam-repseq.fa", sequence)
    write("extension/final_consensuses/fam1_rep.fa", sequence)
    write("extension/images_and_alignments/fam1.fa", sequence)
    for index in range(30):
        write(f"extension/extract_align/fragment-{index}.fa", "AAAA" * 1000)
    summary = write("extension/extension_summary.tsv", "te_id\tstatus\nfam1\tok\n")
    write("extension/extension.manifest.json", json.dumps({"summary_tsv": str(summary), "n_ok": 1, "n_failed": 0}))
    original = write("reclassify/Test_extended_rep.fa", sequence)
    classified = write("reclassify/Test_extended_rep.fa.classified", sequence)
    write("reclassify/repeatclassifier.checkpoint.json", json.dumps(
        reclassify.checkpoint_data(original, classified, "RepeatClassifier")))
    curated = write("curation/Final.RepeatModeler.Lib.fa", sequence)
    write("curation/curation.manifest.json", json.dumps({"final_library": str(curated)}))
    structural = write("classification/Test.structure_input.fa", sequence)
    final_library = write("classification/Test.classified.fa", ">fam1#RC/Helitron\nACGT\n")
    outputs = {}
    for key in ("structure_evidence", "classifications_tsv", "evidence_tsv"):
        outputs[key] = str(write(f"classification/{key}.tsv", "complete\n"))
    outputs.update(input_library=str(curated), classified_library=str(final_library), input_signature={
        "library_sha256": classification.sha256_file(curated),
        "structure_library_sha256": classification.sha256_file(structural),
        "extensionwork_sha256": classification.selected_extensionwork_sha256(root / "extension/extensionwork"),
    })
    write("classification/classification.manifest.json", json.dumps(outputs))
    for index in range(5):
        write(f"classification/dfam/chunks/chunk_{index}.fa", sequence)
    final = []
    for suffix in (".out", ".tbl", ".masked", ".cat.gz", ".align"):
        final.append(write(f"final_annotation/Test.fa{suffix}", f"final{suffix}\n" * 1000))
    write("final_annotation/final_annotation.manifest.json", json.dumps({"repeatmasker_out": str(final[0])}))
    write_checkpoint(root / "final_annotation/repeatmasker.checkpoint.json",
                     input_signature(genome, final_library,
                                     ["RepeatMasker", "-s", "-xsmall", "-gc", "-gff", "-norna", "-a"]), final)
    for name in ("annotation_refinement.tsv", "filteredRepeats.tsv", "annotation_refinement.stats.tsv"):
        write("Defragmenting_RModAnnotation/Test." + name, "complete\n")
    write("Defragmenting_RModAnnotation/Test.annotation_refinement.manifest.json",
          json.dumps({"refinement_algorithm_version": 1}))
    for label in ("noFiltered", "Filtered"):
        write(f"Test.Defragmenting_RModAnnotation_{label}/Test.highLevelCount.txt", "summary\n")
    return config


def test_completed_run_compact_restore_preserves_checkpoint_inputs(completed_run):
    config = completed_run
    root = config.outdir_path
    before_signature = classification.selected_extensionwork_sha256(root / "extension/extensionwork")
    final_hash = classification.sha256_file(root / "final_annotation/Test.fa.out")
    before_mtime = (root / "final_annotation/Test.fa.align").stat().st_mtime_ns
    validate_completed_run(config)
    with output_lock(root):
        report = compact_run(config, logging.getLogger("test.storage"))
    assert report["after"]["files_and_symlinks"] < report["before"]["files_and_symlinks"]
    assert report["after"]["allocated_bytes"] < report["before"]["allocated_bytes"]
    assert not (root / "extension/extensionwork").exists()
    assert not (root / "discovery/rmodeler_dir/RM_old").exists()
    assert not (root / "final_annotation/Test.fa.align").exists()
    assert (root / "final_annotation/Test.fa.out").is_file()
    assert (root / "classification/Test.classified.fa").is_file()
    assert classification.sha256_file(root / "final_annotation/Test.fa.out") == final_hash
    assert compact_run(config, logging.getLogger("test.storage"))["actions"] == []
    with output_lock(root):
        restored = ArchiveStore(root).restore_all()
    assert "extension/extensionwork" in restored
    assert classification.selected_extensionwork_sha256(root / "extension/extensionwork") == before_signature
    assert (root / "final_annotation/Test.fa.align").stat().st_mtime_ns == before_mtime
    validate_completed_run(config)


def test_legacy_repeatmasker_adoption_survives_compact_restore(completed_run):
    from drayte.pipeline import final_annotation

    config = completed_run
    root = config.outdir_path
    marker = root / "final_annotation/repeatmasker.checkpoint.json"
    marker.unlink()
    with output_lock(root):
        compact_run(config, logging.getLogger("test.storage"))
        assert json.loads(marker.read_text())["legacy_adopted"]
        ArchiveStore(root).restore_all()
    library = root / "classification/Test.classified.fa"
    library.write_text(library.read_text())
    # A missing RepeatMasker executable would fail if restore caused a rerun.
    result = final_annotation.run(config, {"classified_library": str(library)}, logging.getLogger("test.storage"))
    assert Path(result["repeatmasker_out"]).is_file()


def test_incomplete_run_keeps_scratch_and_evidence(completed_run):
    config = completed_run
    (config.outdir_path / "final_annotation/Test.fa.out").unlink()
    with pytest.raises(ValueError, match="missing/empty"):
        compact_run(config, logging.getLogger("test.storage"))
    assert (config.outdir_path / "discovery/rmodeler_dir/RM_old").is_dir()
    assert (config.outdir_path / "extension/extensionwork").is_dir()


def test_changed_classification_input_prevents_compaction(completed_run):
    config = completed_run
    (config.outdir_path / "extension/extensionwork/fam1/repam-ranges.tsv").write_text("changed\n")
    with pytest.raises(ValueError, match="fingerprints"):
        compact_run(config, logging.getLogger("test.storage"))
    assert (config.outdir_path / "discovery/rmodeler_dir/RM_old").is_dir()


def test_changed_repeatmasker_output_prevents_compaction(completed_run):
    config = completed_run
    (config.outdir_path / "final_annotation/Test.fa.align").write_text("changed\n")
    with pytest.raises(ValueError, match="RepeatMasker completion"):
        compact_run(config, logging.getLogger("test.storage"))
    assert (config.outdir_path / "discovery/rmodeler_dir/RM_old").is_dir()


def test_scratch_dependency_via_symlink_or_original_genome_is_protected(completed_run):
    config = completed_run
    target = config.outdir_path / "discovery/rmodeler_dir/RM_old"
    link = config.outdir_path / "reference.fa"
    link.symlink_to(target / "round-1/work.fa")
    with pytest.raises(ValueError, match="retained symlink"):
        assert_safe_deletions(config, [target])
    link.unlink()
    config.genome = str(target / "round-1/work.fa")
    with pytest.raises(ValueError, match="configured input"):
        assert_safe_deletions(config, [target])


def test_stage_scratch_cleanup_requires_success_and_compact_opt_in(completed_run):
    config = completed_run
    config.extra["storage"]["mode"] = "keep"
    assert cleanup_stage_scratch(config, "discovery", logging.getLogger("test")) == []
    config.extra["storage"]["mode"] = "compact"
    table = config.outdir_path / "discovery/rmasker_dir/Test.fa.tbl"
    table.unlink()
    with pytest.raises(ValueError, match="missing/empty"):
        cleanup_stage_scratch(config, "discovery", logging.getLogger("test"))
    assert (config.outdir_path / "discovery/rmodeler_dir/RM_old").is_dir()


def test_lock_prevents_overlapping_pipeline_or_cleanup(completed_run):
    with output_lock(completed_run.outdir_path):
        with pytest.raises(RuntimeError, match="holds the output lock"):
            with output_lock(completed_run.outdir_path):
                pass


def test_pipeline_restores_before_any_stage_and_compacts_after_success(completed_run, tmp_path, monkeypatch):
    import yaml

    config = completed_run
    root = config.outdir_path
    ArchiveStore(root).archive("extension/extensionwork")
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({"genome": config.genome, "outdir": config.outdir,
                                    "species": config.species, **config.extra}))
    monkeypatch.setattr(sys, "argv", ["drayte", "--config", str(path)])
    seen = []

    def pipeline(loaded, filename):
        assert (root / "extension/extensionwork/fam1/repam-ranges.tsv").is_file()
        validate_completed_run(loaded)
        seen.append("pipeline")

    monkeypatch.setattr(runner, "run_pipeline", pipeline)
    monkeypatch.setattr(runner, "setup_logger", lambda **kwargs: logging.getLogger("test.storage"))
    runner.main()
    assert seen == ["pipeline"]
    assert not (root / "extension/extensionwork").exists()
    assert ArchiveStore(root).records()


def test_failed_pipeline_retains_restored_evidence(completed_run, tmp_path, monkeypatch):
    import yaml

    config = completed_run
    root = config.outdir_path
    ArchiveStore(root).archive("extension/extensionwork")
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({"genome": config.genome, "outdir": config.outdir,
                                    "species": config.species, **config.extra}))
    monkeypatch.setattr(sys, "argv", ["drayte", "--config", str(path)])
    def fail(*args):
        raise RuntimeError("simulated computational failure")
    monkeypatch.setattr(runner, "run_pipeline", fail)
    with pytest.raises(RuntimeError, match="computational failure"):
        runner.main()
    assert (root / "extension/extensionwork/fam1/repam-ranges.tsv").is_file()
    assert ArchiveStore(root).records() == {}


def test_compact_family_inspection_avoids_duplicate_archives(completed_run, monkeypatch):
    from drayte.reporting import family_inspection

    config = completed_run
    config.extra["family_inspection"] = {"enabled": True}
    seen = []
    def build(**kwargs):
        seen.append(kwargs["archive"])
        return {"families": 1}
    monkeypatch.setattr(family_inspection, "build_family_inspection", build)
    family_inspection.run(config, {"filtered_tsv": str(config.outdir_path / "Defragmenting_RModAnnotation/Test.filteredRepeats.tsv")},
                          logging.getLogger("test"))
    assert seen == [False]


def test_storage_cli_preview_is_read_only_and_shows_validation(completed_run, tmp_path, capsys):
    import yaml
    from drayte.storage.cli import main

    config = completed_run
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({"genome": config.genome, "outdir": config.outdir,
                                    "species": config.species, **config.extra}))
    assert main(["--config", str(path)]) == 0
    assert '"validation": "completed outputs validated"' in capsys.readouterr().out
    assert not (config.outdir_path / ".storage").exists()


def test_old_run_preserves_unknown_cli_only_input_in_scratch(completed_run):
    config = completed_run
    root = config.outdir_path
    imported = root / "discovery/rmodeler_dir/RM_old/imported.fa"
    imported.write_text(">original-cli-input\nACGT\n")
    original = imported.read_bytes()
    # The old YAML does not record the original --repeatmodeler-library input.
    with output_lock(root):
        report = compact_run(config, logging.getLogger("test.storage"))
        assert not any(action["operation"] == "delete_scratch" for action in report["actions"])
        assert "discovery/rmodeler_dir/RM_old" in ArchiveStore(root).records()
        ArchiveStore(root).restore_all()
    assert imported.read_bytes() == original


@pytest.mark.parametrize("settings", [[], False, "compact", {"mode": "invalid"}])
def test_invalid_storage_settings_rejected(completed_run, settings):
    completed_run.extra["storage"] = settings
    with pytest.raises(ValueError):
        storage_mode(completed_run)


def test_plan_does_not_mutate_completed_output(completed_run):
    root = completed_run.outdir_path
    assert not (root / ".storage").exists()
    assert any(item["path"] == "extension/extensionwork" for item in plan_compaction(completed_run))
    assert not (root / ".storage").exists()


def prepare_early_extension(config):
    root = config.outdir_path
    work = root / "extension/extensionwork/fam1"
    sequence = (root / "extension/final_consensuses/fam1_rep.fa").read_text()
    (work / "fam1_rep.fa").write_text(sequence)
    (work / "rep").write_text(sequence)
    (work / "repseq.unextended").write_text(">copy1\nACGT\n" * 12)
    (work / "MSA-extended_with_rmod_cons.fa").write_text(sequence)
    (work / "img.png").write_bytes(b"image")
    (work / "scratch.fa").write_text("ACGT" * 50000)
    return read_extension_manifest(config), work


def read_extension_manifest(config):
    return json.loads((config.outdir_path / "extension/extension.manifest.json").read_text())


def test_early_extension_compacts_before_final_stage_and_preserves_consumers(completed_run):
    from drayte.extension.consensus import postprocess_extension_outputs

    config = completed_run
    result, work = prepare_early_extension(config)
    final_output = config.outdir_path / "final_annotation/Test.fa.out"
    final_output.unlink()  # Whole pipeline is unfinished; extension is complete.
    before = classification.selected_extensionwork_sha256(work.parent)
    expected = postprocess_extension_outputs("fam1", work)[3:]
    with output_lock(config.outdir_path):
        actions = compact_completed_extension(config, result, logging.getLogger("test"))
        assert actions
        assert not (work / "scratch.fa").exists()
        assert (work / "fam1_rep.fa").is_file()
        assert (work / "repam-ranges.tsv").is_file()
        assert classification.selected_extensionwork_sha256(work.parent) == before
        assert postprocess_extension_outputs("fam1", work)[3:] == expected == (12, 4)
        assert compact_completed_extension(config, result, logging.getLogger("test")) == []
        ArchiveStore(config.outdir_path).restore_all()
    assert (work / "scratch.fa").is_file()


def test_final_compaction_accepts_early_extension_archive(completed_run):
    config = completed_run
    result, work = prepare_early_extension(config)
    with output_lock(config.outdir_path):
        compact_completed_extension(config, result, logging.getLogger("test"))
        validate_completed_run(config)
        report = compact_run(config, logging.getLogger("test"))
        assert report["actions"]
        assert "extension/extensionwork" in report["archives"]
        ArchiveStore(config.outdir_path).restore_all()
        validate_completed_run(config)
    assert (work / "scratch.fa").is_file()


def test_incomplete_extension_keeps_every_file(completed_run):
    config = completed_run
    result, work = prepare_early_extension(config)
    result["n_failed"] = 1
    assert compact_completed_extension(config, result, logging.getLogger("test")) == []
    assert (work / "scratch.fa").is_file()
    assert ArchiveStore(config.outdir_path).records() == {}


def test_early_extension_preserves_external_consumers(completed_run):
    config = completed_run
    result, work = prepare_early_extension(config)
    linked = work / "extra_evidence.tsv"
    linked.write_text("evidence\n")
    (config.outdir_path / "linked.tsv").symlink_to(linked)
    supplied = work / "custom_input.fa"
    supplied.write_text(">f\nACGT\n")
    config.extra["other_input"] = str(supplied)
    compact_completed_extension(config, result, logging.getLogger("test"))
    assert linked.is_file() and supplied.is_file()
    assert not (work / "scratch.fa").exists()


def test_early_extension_keeps_repam_from_unlisted_legacy_families(completed_run):
    config = completed_run
    result, work = prepare_early_extension(config)
    legacy = work.parent / "legacy-family"
    legacy.mkdir()
    (legacy / "repam-ranges.tsv").write_text("1\t4\n")
    (legacy / "repam-repseq.fa").write_text(">old\nAAAA\n")
    before = classification.selected_extensionwork_sha256(work.parent)
    compact_completed_extension(config, result, logging.getLogger("test"))
    assert classification.selected_extensionwork_sha256(work.parent) == before
    assert (legacy / "repam-repseq.fa").is_file()


def test_explicit_default_repam_directory_allows_early_compaction(completed_run):
    config = completed_run
    result, work = prepare_early_extension(config)
    config.extra["classification"] = {"repam_extensionwork_dir": str(work.parent)}
    compact_completed_extension(config, result, logging.getLogger("test"))
    assert not (work / "scratch.fa").exists()
    assert config.extra["classification"]["repam_extensionwork_dir"] == str(work.parent)
