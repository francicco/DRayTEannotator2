import csv

from drayte.pipeline.extension import validate_extension_outputs


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
