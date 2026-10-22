import csv
from types import SimpleNamespace

from drayte.pipeline.extension import (
    candidate_matches_cached_extension_query,
    extension_signatures_compatible,
    fasta_libraries_equivalent,
    trust_complete_legacy_checkpoint,
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
    write_summary(
        summary,
        [
            {
                "te_id": family_id,
                "status": "ok",
            }
        ],
    )

    valid, reason = validate_extension_outputs(
        {
            "summary_tsv": str(summary),
            "n_ok": 1,
            "n_failed": 0,
        },
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
        {
            "summary_tsv": str(summary),
            "n_ok": 1,
            "n_failed": 0,
        },
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
        {
            "summary_tsv": str(summary),
            "n_ok": 0,
            "n_failed": 1,
        },
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


def test_fasta_libraries_equivalent_ignores_order_wrapping_and_descriptions(tmp_path):
    current = tmp_path / "current.fa"
    cached = tmp_path / "cached.fa"

    current.write_text(
        ">family-1#Unknown current description\nACGTACGT\n"
        ">family-2#LINE/L1\nTTTTAAAA\n"
    )
    cached.write_text(
        ">family-2#LINE/L1 cached description\nTTTT\nAAAA\n"
        ">family-1#Unknown\nACGT\nACGT\n"
    )

    equivalent, reason = fasta_libraries_equivalent(current, cached)

    assert equivalent
    assert reason == "same 2 identifiers and sequences"


def test_fasta_libraries_equivalent_rejects_changed_sequence(tmp_path):
    current = tmp_path / "current.fa"
    cached = tmp_path / "cached.fa"
    current.write_text(">family-1#Unknown\nACGT\n")
    cached.write_text(">family-1#Unknown\nACGA\n")

    equivalent, reason = fasta_libraries_equivalent(current, cached)

    assert not equivalent
    assert reason == "1 family sequences differ"


def test_fasta_libraries_equivalent_rejects_changed_identifiers(tmp_path):
    current = tmp_path / "current.fa"
    cached = tmp_path / "cached.fa"
    current.write_text(">family-new#Unknown\nACGT\n")
    cached.write_text(">family-old#Unknown\nACGT\n")

    equivalent, reason = fasta_libraries_equivalent(current, cached)

    assert not equivalent
    assert "cached identifiers are missing" in reason
    assert "current identifiers are absent" in reason


def test_candidate_matches_legacy_filtered_and_encoded_extension_query(tmp_path):
    candidate = tmp_path / "external-library.fa"
    cached = tmp_path / "cached-extension-query.fa"
    candidate.write_text(
        ">rnd-1_family-1#LINE/L1\n" + "ACGT" * 30 + "\n"
        ">rnd-1_family-2#Unknown\n" + "TGCA" * 28 + "\n"
        ">rnd-1_family-short#Unknown\n" + "A" * 99 + "\n"
    )
    cached.write_text(
        ">Mmim-rnd-1_family-2__Unknown\n" + "TGCA" * 28 + "\n"
        ">Mmim-rnd-1_family-1__LINE___L1\n" + "ACGT" * 30 + "\n"
    )

    equivalent, reason = candidate_matches_cached_extension_query(
        candidate_library=candidate,
        cached_query=cached,
        species="Mmim",
    )

    assert equivalent
    assert "same 2 extension-query identifiers and sequences" in reason


def test_candidate_rejects_changed_legacy_extension_sequence(tmp_path):
    candidate = tmp_path / "external-library.fa"
    cached = tmp_path / "cached-extension-query.fa"
    candidate.write_text(">rnd-1_family-1#Unknown\n" + "A" * 100 + "\n")
    cached.write_text(">Mmim-rnd-1_family-1__Unknown\n" + "C" * 100 + "\n")

    equivalent, reason = candidate_matches_cached_extension_query(
        candidate_library=candidate,
        cached_query=cached,
        species="Mmim",
    )

    assert not equivalent
    assert reason == "1 family sequences differ"


def test_trusted_legacy_checkpoint_requires_complete_outputs_and_pins_query(tmp_path):
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
    config = SimpleNamespace(extra={"classification": {}})
    result = {
        "summary_tsv": str(summary),
        "n_ok": 1,
        "n_failed": 0,
    }

    trusted, reason = trust_complete_legacy_checkpoint(
        result,
        extensionwork=extensionwork,
        final_consensuses=final_consensuses,
        cached_query=cached_query,
        config=config,
    )

    assert trusted, reason
    assert config.extra["classification"]["structure_base_fasta"] == str(cached_query)
    assert result["resume_validation"]["mode"] == "trusted_legacy_extension"
    assert "input_signature" not in result


def test_trusted_legacy_checkpoint_rejects_incomplete_outputs(tmp_path):
    extensionwork = tmp_path / "extensionwork"
    final_consensuses = tmp_path / "final_consensuses"
    extensionwork.mkdir()
    final_consensuses.mkdir()
    summary = tmp_path / "extension_summary.tsv"
    write_summary(summary, [{"te_id": "missing-family", "status": "ok"}])
    cached_query = tmp_path / "cached.fa"
    cached_query.write_text(">missing-family\nACGT\n")
    config = SimpleNamespace(extra={})

    trusted, reason = trust_complete_legacy_checkpoint(
        {
            "summary_tsv": str(summary),
            "n_ok": 1,
            "n_failed": 0,
        },
        extensionwork=extensionwork,
        final_consensuses=final_consensuses,
        cached_query=cached_query,
        config=config,
    )

    assert not trusted
    assert "missing extension work directory" in reason
    assert "classification" not in config.extra
