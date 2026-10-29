import hashlib
import json
import logging
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from drayte.classification.structure_inputs import resolve_structure_input_fasta
from drayte.pipeline import classification as classification_stage
from drayte.pipeline.extension import prepare_complete_extension_checkpoint
from drayte.pipeline.classification import (
    classification_input_signature,
    classification_signatures_compatible,
    classified_library_matches_input,
)


def test_classified_library_matches_input_ignores_classification_suffix(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">family-1#Unknown\nACGT\n")

    classified = tmp_path / "classified.fa"
    classified.write_text(">family-1#DNA/TcMar\nACGT\n")

    assert classified_library_matches_input(source, classified)


def test_classified_library_does_not_match_changed_sequence(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">family-1#Unknown\nACGT\n")

    classified = tmp_path / "classified.fa"
    classified.write_text(">family-1#DNA/TcMar\nAAAA\n")

    assert not classified_library_matches_input(source, classified)


def make_signature(tmp_path, classification_cfg):
    library = tmp_path / "library.fa"
    structure_library = tmp_path / "structure.fa"
    extensionwork = tmp_path / "extensionwork"
    extensionwork.mkdir(exist_ok=True)
    library.write_text(">family-1#Unknown\nACGT\n")
    structure_library.write_text(">family-1#Unknown\nACGT\n")
    return classification_input_signature(
        library=library,
        structure_library=structure_library,
        extensionwork_dir=extensionwork,
        classification_cfg=classification_cfg,
        stage_name="classification",
        final_mode=True,
    )


def test_execution_only_settings_do_not_invalidate_classification(tmp_path):
    cached = make_signature(
        tmp_path,
        {
            "min_orf_nt": 300,
            "chunks": 40,
            "jobs": 40,
            "cpu": 1,
        },
    )
    current = make_signature(
        tmp_path,
        {
            "min_orf_nt": 300,
            "chunks": 10,
            "jobs": 10,
            "cpu": 4,
        },
    )

    compatible, legacy_migration = classification_signatures_compatible(
        cached,
        current,
    )

    assert compatible
    assert not legacy_migration
    assert cached["execution_config"] != current["execution_config"]


def test_result_affecting_setting_invalidates_classification(tmp_path):
    cached = make_signature(tmp_path, {"min_orf_nt": 300, "chunks": 40})
    current = make_signature(tmp_path, {"min_orf_nt": 600, "chunks": 10})

    compatible, legacy_migration = classification_signatures_compatible(
        cached,
        current,
    )

    assert not compatible
    assert not legacy_migration


def test_legacy_signature_with_same_biological_inputs_can_be_migrated(tmp_path):
    current = make_signature(tmp_path, {"min_orf_nt": 300, "chunks": 10})
    legacy = {
        "schema_version": 1,
        **{key: current[key] for key in (
            "stage_name",
            "final_mode",
            "library_sha256",
            "structure_library_sha256",
            "extensionwork_sha256",
        )},
        "classification_config_sha256": "legacy-full-config-hash",
    }

    compatible, legacy_migration = classification_signatures_compatible(
        legacy,
        current,
    )

    assert compatible
    assert legacy_migration


def legacy_v2_signature(signature, classification_cfg):
    """Reproduce the previous signature format independently of current code."""
    result_config = {
        key: value for key, value in classification_cfg.items()
        if key not in {"chunks", "jobs", "cpu"}
    }
    return {
        "schema_version": 2,
        **{key: signature[key] for key in (
            "stage_name", "final_mode", "library_sha256",
            "structure_library_sha256", "extensionwork_sha256",
        )},
        "result_config_sha256": hashlib.sha256(
            json.dumps(result_config, sort_keys=True, default=str).encode("utf-8")
        ).hexdigest(),
        "execution_config": {key: classification_cfg[key] for key in ("chunks", "jobs", "cpu") if key in classification_cfg},
    }


@pytest.fixture
def extension_checkpoint(tmp_path):
    root = tmp_path / "run"
    discovery = root / "discovery" / "rmodeler_dir"
    discovery.mkdir(parents=True)
    base = discovery / "T-families.mod.fa"
    base.write_text(">T-family-1#Unknown\nAAAA\n>T-family-2#Unknown\nCCCC\n")
    extension = root / "extension"
    family_dir = extension / "extensionwork" / "T-family-1#Unknown"
    family_dir.mkdir(parents=True)
    (family_dir / "repam-ranges.tsv").write_text("chr1\t1\t8\n")
    (family_dir / "repam-repseq.fa").write_text(">copy-1\nACGTACGT\n")
    final_consensuses = extension / "final_consensuses"
    final_consensuses.mkdir()
    (final_consensuses / "T-family-1#Unknown_rep.fa").write_text(
        ">T-family-1#Unknown\nACGTACGT\n"
    )
    cached_query = extension / "blastfiles" / base.name
    cached_query.parent.mkdir()
    cached_query.write_text(base.read_text())
    summary = extension / "extension_summary.tsv"
    summary.write_text("te_id\tstatus\nT-family-1#Unknown\tok\n")
    extension_result = {"summary_tsv": str(summary), "n_ok": 1, "n_failed": 0}
    library = root / "curated.fa"
    library.write_text(">T-family-1#Unknown\nACGTACGT\n")
    config = SimpleNamespace(
        species="T", outdir_path=root, threads=2,
        extra={"classification": {
            "min_orf_nt": 300, "pfam_db": "/db/Pfam-A.hmm",
            "chunks": 40, "jobs": 4, "cpu": 1,
            # Force resolution to read the selected base on both calls.
            "rebuild_structure_input_fasta": True,
        }},
    )
    outdir = root / "classification"
    outdir.mkdir()

    def signature():
        structure_library = resolve_structure_input_fasta(
            config=config,
            classification_cfg=config.extra["classification"],
            outdir=outdir,
            fallback_fasta=library,
        )
        return classification_input_signature(
            library=library,
            structure_library=structure_library,
            extensionwork_dir=extension / "extensionwork",
            classification_cfg=config.extra["classification"],
            stage_name="classification", final_mode=True,
        )

    def pin_checkpoint():
        valid, reason = prepare_complete_extension_checkpoint(
            extension_result,
            extensionwork=extension / "extensionwork",
            final_consensuses=final_consensuses,
            cached_query=cached_query,
            config=config,
        )
        assert valid, reason

    return SimpleNamespace(
        config=config, library=library, base=base, cached_query=cached_query,
        family_dir=family_dir, outdir=outdir, signature=signature,
        pin_checkpoint=pin_checkpoint,
        structure_library=outdir / "T.structure_input.fa",
    )


def test_resuming_extension_does_not_change_classification_signature(extension_checkpoint):
    case = extension_checkpoint
    fresh = case.signature()
    fresh_structure = case.structure_library.read_bytes()
    assert "structure_base_fasta" not in case.config.extra["classification"]

    case.pin_checkpoint()
    resumed = case.signature()

    assert case.config.extra["classification"]["structure_base_fasta"] == str(case.cached_query)
    assert case.structure_library.read_bytes() == fresh_structure
    assert fresh["result_config_sha256"] == resumed["result_config_sha256"]
    assert fresh["schema_version"] == resumed["schema_version"] == 3
    assert classification_signatures_compatible(fresh, resumed) == (True, False)


@pytest.mark.parametrize("cached_with_path", [False, True])
def test_v2_signature_migrates_with_or_without_injected_structure_path(
    extension_checkpoint, cached_with_path
):
    case = extension_checkpoint
    if cached_with_path:
        case.pin_checkpoint()
    cached = legacy_v2_signature(case.signature(), case.config.extra["classification"])
    case.pin_checkpoint()
    case.config.extra["classification"].update(chunks=10, jobs=2, cpu=4)
    current = case.signature()

    assert classification_signatures_compatible(cached, current) == (True, True)


@pytest.mark.parametrize("version", [2, 3])
@pytest.mark.parametrize("key, changed", [
    ("min_orf_nt", 600),
    ("pfam_db", "/db/other-Pfam.hmm"),
    ("dfam_db", "/db/Dfam.hmm"),
    ("forward_only", True),
    ("mmseqs_min_identity", 0.95),
    ("unrelated_custom_setting", "must-still-affect-the-hash"),
])
def test_structure_path_normalization_does_not_ignore_other_settings(
    extension_checkpoint, version, key, changed
):
    case = extension_checkpoint
    cached = case.signature()
    if version == 2:
        cached = legacy_v2_signature(cached, case.config.extra["classification"])
    case.pin_checkpoint()
    case.config.extra["classification"][key] = changed

    assert classification_signatures_compatible(cached, case.signature()) == (False, False)


@pytest.mark.parametrize("version", [2, 3])
@pytest.mark.parametrize("changed_input", ["curated", "structure", "extensionwork"])
def test_changed_biological_input_rejects_path_normalized_checkpoint(
    extension_checkpoint, version, changed_input
):
    case = extension_checkpoint
    cached = case.signature()
    if version == 2:
        cached = legacy_v2_signature(cached, case.config.extra["classification"])
    case.pin_checkpoint()
    if changed_input == "curated":
        case.library.write_text(">T-family-1#Unknown\nACGTAAAA\n")
    elif changed_input == "structure":
        case.cached_query.write_text(
            ">T-family-1#Unknown\nAAAA\n>T-family-2#Unknown\nGGGG\n"
        )
    else:
        (case.family_dir / "repam-repseq.fa").write_text(">copy-1\nTGTGTGTG\n")

    assert classification_signatures_compatible(cached, case.signature()) == (False, False)


def test_unknown_signature_schema_is_rejected(extension_checkpoint):
    current = extension_checkpoint.signature()
    assert classification_signatures_compatible(
        {**current, "schema_version": 999}, current
    ) == (False, False)


def test_v2_migration_requires_a_matching_historical_result_config_hash(extension_checkpoint):
    current = extension_checkpoint.signature()
    cached = legacy_v2_signature(current, extension_checkpoint.config.extra["classification"])
    cached["result_config_sha256"] = "unrecognized-configuration-hash"
    assert classification_signatures_compatible(cached, current) == (False, False)


def write_classification_checkpoint(case, signature):
    result = {"stage": "classification", "input_signature": signature}
    for key, filename in (
        ("structure_evidence", "T.structure_summary.tsv"),
        ("classifications_tsv", "T.classifications.tsv"),
        ("evidence_tsv", "T.evidence.tsv"),
    ):
        path = case.outdir / filename
        path.write_text("complete cached output\n")
        result[key] = str(path)
    classified = case.outdir / "T.classified.fa"
    classified.write_text(case.structure_library.read_text().replace("#Unknown", "#DNA/TcMar"))
    result["classified_library"] = str(classified)
    manifest = case.outdir / "classification.manifest.json"
    manifest.write_text(json.dumps(result))
    return result, manifest


@pytest.mark.parametrize("version", ["v3", "v2-no-path", "v2-with-path"])
def test_classification_run_reuses_checkpoint_after_extension_resume(
    extension_checkpoint, monkeypatch, version
):
    case = extension_checkpoint
    if version == "v2-with-path":
        case.pin_checkpoint()
    signature = case.signature()
    if version != "v3":
        signature = legacy_v2_signature(signature, case.config.extra["classification"])
    cached, manifest = write_classification_checkpoint(case, signature)
    case.pin_checkpoint()

    def unexpected(*args, **kwargs):
        pytest.fail("validated classification checkpoint must skip evidence generation")

    monkeypatch.setattr(classification_stage, "detect_complexity_from_fasta", unexpected)
    monkeypatch.setattr(classification_stage, "run_domain_annotation", unexpected)
    result = classification_stage.run(
        case.config, {"final_library": str(case.library)},
        logging.getLogger("test.classification.resume"), final_mode=True,
    )

    assert result["classified_library"] == cached["classified_library"]
    assert result["input_signature"]["schema_version"] == 3
    assert json.loads(manifest.read_text())["input_signature"] == result["input_signature"]


def test_v2_migration_rejects_classified_library_with_changed_sequence(
    extension_checkpoint,
):
    case = extension_checkpoint
    signature = legacy_v2_signature(case.signature(), case.config.extra["classification"])
    result, _ = write_classification_checkpoint(case, signature)
    Path(result["classified_library"]).write_text(
        ">T-family-1#DNA/TcMar\nTTTTTTTT\n>T-family-2#DNA/TcMar\nCCCC\n"
    )
    case.pin_checkpoint()
    with pytest.raises(RuntimeError, match="cannot be safely migrated"):
        classification_stage.run(
            case.config, {"final_library": str(case.library)},
            logging.getLogger("test.classification.resume.mismatch"), final_mode=True,
        )


def test_v2_migration_validates_all_required_outputs_before_skipping(
    extension_checkpoint, monkeypatch
):
    case = extension_checkpoint
    signature = legacy_v2_signature(case.signature(), case.config.extra["classification"])
    result, _ = write_classification_checkpoint(case, signature)
    Path(result["evidence_tsv"]).unlink()
    case.pin_checkpoint()

    class ResumeEvidenceGeneration(Exception):
        pass

    def resumed(*args, **kwargs):
        raise ResumeEvidenceGeneration

    monkeypatch.setattr(classification_stage, "detect_complexity_from_fasta", resumed)
    with pytest.raises(ResumeEvidenceGeneration):
        classification_stage.run(
            case.config, {"final_library": str(case.library)},
            logging.getLogger("test.classification.resume.incomplete"), final_mode=True,
        )


@pytest.mark.parametrize("key, changed", [
    ("min_orf_nt", 600),
    ("pfam_db", "/db/changed-Pfam.hmm"),
])
def test_signature_mismatch_diagnostic_names_changed_setting(
    extension_checkpoint, key, changed
):
    case = extension_checkpoint
    cached = case.signature()
    case.pin_checkpoint()
    case.config.extra["classification"][key] = changed
    current = case.signature()

    assert classification_stage.classification_signature_mismatch_details(
        cached, current
    ) == f"classification.{key}"


@pytest.mark.parametrize("version", [2, 3])
@pytest.mark.parametrize("changed_sequence", [False, True])
def test_checkpoint_checks_older_alternate_structure_source_content(
    extension_checkpoint, monkeypatch, version, changed_sequence
):
    case = extension_checkpoint
    # Exercise the normal stage behavior without a user-requested rebuild.
    case.config.extra["classification"].pop("rebuild_structure_input_fasta")
    signature = case.signature()
    if version == 2:
        signature = legacy_v2_signature(signature, case.config.extra["classification"])
    cached, _ = write_classification_checkpoint(case, signature)
    old_structure = case.structure_library.read_bytes()

    alternate_base = case.config.outdir_path / "older-alternate-base.fa"
    alternate_text = case.base.read_text()
    if changed_sequence:
        # Family 2 has no extension override, so the selected base determines
        # its biological sequence even though it is absent from curated.fa.
        alternate_text = alternate_text.replace("CCCC", "GGGG")
    alternate_base.write_text(alternate_text)
    older_ns = case.structure_library.stat().st_mtime_ns - 10_000_000_000
    os.utime(alternate_base, ns=(older_ns, older_ns))
    assert alternate_base.stat().st_mtime_ns < case.structure_library.stat().st_mtime_ns
    case.config.extra["classification"]["structure_base_fasta"] = str(alternate_base)

    def unexpected(*args, **kwargs):
        pytest.fail("source validation must happen before any external evidence work")

    monkeypatch.setattr(classification_stage, "detect_complexity_from_fasta", unexpected)
    monkeypatch.setattr(classification_stage, "run_domain_annotation", unexpected)

    def run():
        return classification_stage.run(
            case.config, {"final_library": str(case.library)},
            logging.getLogger("test.classification.older-structure-source"), final_mode=True,
        )

    if changed_sequence:
        with pytest.raises(RuntimeError, match="structure_library_sha256"):
            run()
        assert case.structure_library.read_bytes() != old_structure
        assert "GGGG" in case.structure_library.read_text()
    else:
        result = run()
        assert result["classified_library"] == cached["classified_library"]
        assert case.structure_library.read_bytes() == old_structure
        assert result["input_signature"]["schema_version"] == 3
