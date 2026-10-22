import csv

from drayte.pipeline.extension import (
    extension_signatures_compatible,
    fasta_libraries_equivalent,
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
