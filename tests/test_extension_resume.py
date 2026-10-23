import csv
import json
import logging
from types import SimpleNamespace

from drayte.pipeline import extension as extension_stage
from drayte.pipeline.extension import (
    extension_signatures_compatible,
    prepare_complete_extension_checkpoint,
    validate_extension_outputs,
)


def write_summary(path, rows):
    fieldnames = [
        "te_id",
        "cat_file",
        "status",
        "category",
        "rep_fa",
        "msa_fa",
        "png_file",
        "hit_count",
        "consensus_length",
        "error_message",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_validate_extension_outputs_accepts_complete_successful_family(tmp_path):
    extensionwork = tmp_path / "extensionwork"
    final_consensuses = tmp_path / "final_consensuses"
    family_id = "Pmex-rnd-1_family-0#Unknown"
    (extensionwork / family_id).mkdir(parents=True)
    final_consensuses.mkdir()
    (final_consensuses / f"{family_id}_rep.fa").write_text(">family\nACGT\n")

    summary = tmp_path / "extension_summary.tsv"
    write_summary(summary, [{"te_id": family_id, "status": "ok"}])

    valid, reason = validate_extension_outputs(
        {"summary_tsv": str(summary), "n_ok": 1, "n_failed": 0},
        extensionwork=extensionwork,
        final_consensuses=final_consensuses,
    )

    assert valid, reason


def test_validate_extension_outputs_rejects_missing_final_consensus(tmp_path):
    extensionwork = tmp_path / "extensionwork"
    final_consensuses = tmp_path / "final_consensuses"
    family_id = "Pmex-rnd-1_family-0#Unknown"
    (extensionwork / family_id).mkdir(parents=True)
    final_consensuses.mkdir()

    summary = tmp_path / "extension_summary.tsv"
    write_summary(summary, [{"te_id": family_id, "status": "ok"}])

    valid, reason = validate_extension_outputs(
        {"summary_tsv": str(summary), "n_ok": 1, "n_failed": 0},
        extensionwork=extensionwork,
        final_consensuses=final_consensuses,
    )

    assert not valid
    assert "final consensus" in reason


def test_validate_extension_outputs_retries_failed_families(tmp_path):
    extensionwork = tmp_path / "extensionwork"
    final_consensuses = tmp_path / "final_consensuses"
    extensionwork.mkdir()
    final_consensuses.mkdir()

    summary = tmp_path / "extension_summary.tsv"
    write_summary(
        summary,
        [{"te_id": "Pmex-rnd-1_family-0#Unknown", "status": "failed_extension"}],
    )

    valid, reason = validate_extension_outputs(
        {"summary_tsv": str(summary), "n_ok": 0, "n_failed": 1},
        extensionwork=extensionwork,
        final_consensuses=final_consensuses,
    )

    assert not valid
    assert "1 failed families" in reason


def test_version_one_signature_is_compatible_with_same_core_inputs():
    cached = {
        "schema_version": 1,
        "species": "Pmex",
        "query_fasta_sha256": "query",
        "genome_fasta_sha256": "genome",
        "extension_max_hits": 50,
        "extension_left_buffer": 100,
        "extension_right_buffer": 100,
        "repeatmodeler_extend_script": "/tools/davidExtendConsRAM.pl",
    }
    current = {
        **cached,
        "schema_version": 2,
        "repeatmodeler_extend_script_sha256": "script",
        "ucsc_tools_dir": "/tools/ucsc",
        "repeatmodeler_dir": "/tools/repeatmodeler",
        "repeat_after_me_dir": "/tools/repeat-after-me",
        "phrap_dir": "/tools/phrap",
    }

    assert extension_signatures_compatible(cached, current)


def test_version_two_signature_detects_changed_tool_directory():
    cached = {
        "schema_version": 2,
        "species": "Pmex",
        "query_fasta_sha256": "query",
        "genome_fasta_sha256": "genome",
        "ucsc_tools_dir": "/tools/ucsc-v1",
    }
    current = {**cached, "ucsc_tools_dir": "/tools/ucsc-v2"}

    assert not extension_signatures_compatible(cached, current)


def test_complete_checkpoint_requires_outputs_and_cached_query(tmp_path):
    extensionwork = tmp_path / "extensionwork"
    final_consensuses = tmp_path / "final_consensuses"
    family_id = "Mmim-rnd-1_family-1__Unknown"
    (extensionwork / family_id).mkdir(parents=True)
    final_consensuses.mkdir()
    (final_consensuses / f"{family_id}_rep.fa").write_text(">family\nACGT\n")

    summary = tmp_path / "extension_summary.tsv"
    write_summary(summary, [{"te_id": family_id, "status": "ok"}])
    cached_query = tmp_path / "Mmim-families.mod.fa"
    cached_query.write_text(f">{family_id}\nACGT\n")
    result = {
        "summary_tsv": str(summary),
        "n_ok": 1,
        "n_failed": 0,
    }

    complete, reason = prepare_complete_extension_checkpoint(
        result,
        extensionwork=extensionwork,
        final_consensuses=final_consensuses,
        cached_query=cached_query,
    )

    assert complete, reason
    assert result["resume_validation"]["mode"] == "validated_extension_checkpoint"
    assert "input_signature" not in result


def test_run_automatically_reuses_complete_legacy_checkpoint(tmp_path):
    outdir = tmp_path / "run"
    assemblies_dir = outdir / "discovery" / "assemblies_dir"
    rmodeler_dir = outdir / "discovery" / "rmodeler_dir"
    assemblies_dir.mkdir(parents=True)
    rmodeler_dir.mkdir(parents=True)

    genome = assemblies_dir / "Mmim.fa"
    genome.write_text(">chr1\n" + "ACGT" * 30 + "\n")
    (assemblies_dir / "Mmim.2bit").write_text("cached two-bit\n")

    candidate = rmodeler_dir / "consensi.fa.classified"
    candidate.write_text(">rnd-current#Unknown\n" + "A" * 120 + "\n")
    current_query = rmodeler_dir / "Mmim-families.mod.fa"
    current_query.write_text(">Mmim-rnd-current#Unknown\n" + "A" * 120 + "\n")

    extension_dir = outdir / "extension"
    extensionwork = extension_dir / "extensionwork"
    final_consensuses = extension_dir / "final_consensuses"
    blastfiles = extension_dir / "blastfiles"
    family_id = "Mmim-rnd-cached__Unknown"
    (extensionwork / family_id).mkdir(parents=True)
    final_consensuses.mkdir(parents=True)
    blastfiles.mkdir(parents=True)
    (final_consensuses / f"{family_id}_rep.fa").write_text(">family\nACGT\n")
    cached_query = blastfiles / "Mmim-families.mod.fa"
    cached_query.write_text(f">{family_id}\n" + "C" * 120 + "\n")

    summary = extension_dir / "extension_summary.tsv"
    write_summary(summary, [{"te_id": family_id, "status": "ok"}])
    manifest = extension_dir / "extension.manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "stage": "extension",
                "outdir": str(extension_dir),
                "extended_library_dir": str(final_consensuses),
                "summary_tsv": str(summary),
                "n_ok": 1,
                "n_failed": 0,
            }
        )
    )

    extend_script = tmp_path / "davidExtendConsRAM.pl"
    extend_script.write_text("#!/usr/bin/env perl\n")
    config = SimpleNamespace(
        outdir_path=outdir,
        species="Mmim",
        genome_path=tmp_path / "source-genome.fa",
        threads=1,
        extra={
            "repeatmodeler_extend_script": str(extend_script),
            "ucsc_tools_dir": str(tmp_path / "kent"),
            "repeatmodeler_dir": str(tmp_path / "RepeatModeler"),
            "repeat_after_me_dir": str(tmp_path / "RepeatAfterMe"),
            "phrap_dir": str(tmp_path / "phrap"),
        },
    )

    result = extension_stage.run(
        config,
        {"raw_library": str(candidate)},
        logging.getLogger("test.extension.auto-resume"),
    )

    assert result["resume_validation"]["input_signature_status"] == "legacy_manifest"
