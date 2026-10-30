import json
import logging
from pathlib import Path

import drayte.pipeline.reclassify as reclassify_module
from drayte.pipeline.reclassify import (
    checkpoint_data,
    classified_output_has_complete_family_set,
    classified_output_matches_input,
    compare_classified_output,
    run_repeatclassifier,
)


def test_classified_output_can_change_class_without_changing_family_sequence(tmp_path):
    merged = tmp_path / "merged.fa"
    merged.write_text(
        ">Pmex-rnd-1_family-0#Unknown\nACGT\n"
        ">Pmex-rnd-1_family-1#Unknown\nTTAA\n"
    )

    classified = tmp_path / "merged.fa.classified"
    classified.write_text(
        ">Pmex-rnd-1_family-0#DNA/TcMar\nACGT\n"
        ">Pmex-rnd-1_family-1#LINE/L1\nTTAA\n"
    )

    assert classified_output_matches_input(merged, classified)


def test_classified_output_is_rejected_when_a_sequence_changed(tmp_path):
    merged = tmp_path / "merged.fa"
    merged.write_text(">Pmex-rnd-1_family-0#Unknown\nACGT\n")

    classified = tmp_path / "merged.fa.classified"
    classified.write_text(">Pmex-rnd-1_family-0#DNA/TcMar\nAAAA\n")

    assert not classified_output_matches_input(merged, classified)


def test_compare_classified_output_reports_mismatch_types(tmp_path):
    merged = tmp_path / "merged.fa"
    merged.write_text(
        ">family-1#Unknown\nACGT\n"
        ">family-2#Unknown\nTTAA\n"
    )
    classified = tmp_path / "merged.fa.classified"
    classified.write_text(
        ">family-1#DNA\nAAAA\n"
        ">family-3#LINE\nCCCC\n"
    )

    comparison = compare_classified_output(merged, classified)

    assert comparison == {
        "input_families": 2,
        "output_families": 2,
        "missing_families": 1,
        "unexpected_families": 1,
        "changed_lengths": 0,
        "changed_sequences": 1,
    }


def write_matching_fastas(tmp_path: Path) -> tuple[Path, Path]:
    merged = tmp_path / "merged.fa"
    merged.write_text(">family-1#Unknown\nACGT\n")
    classified = tmp_path / "merged.fa.classified"
    classified.write_text(">family-1#DNA/TcMar\nACGT\n")
    return merged, classified


def test_completed_checkpoint_reuses_repeatclassifier_output(
    tmp_path,
    monkeypatch,
):
    merged, classified = write_matching_fastas(tmp_path)
    checkpoint = tmp_path / "repeatclassifier.checkpoint.json"
    checkpoint.write_text(
        json.dumps(checkpoint_data(merged, classified, "RepeatClassifier"))
    )

    def unexpected_popen(*args, **kwargs):
        raise AssertionError("RepeatClassifier should have been skipped")

    monkeypatch.setattr(reclassify_module.subprocess, "Popen", unexpected_popen)

    result = run_repeatclassifier(
        merged,
        tmp_path,
        "RepeatClassifier",
        logging.getLogger("test"),
    )

    assert result == classified


def test_completed_checkpoint_reuses_transformed_repeatclassifier_output(
    tmp_path,
    monkeypatch,
):
    merged = tmp_path / "merged.fa"
    merged.write_text(">family-1#Unknown\nACGT\n")
    classified = tmp_path / "merged.fa.classified"
    classified.write_text(">family-1#DNA/TcMar\nAAAA\n")
    checkpoint = tmp_path / "repeatclassifier.checkpoint.json"
    checkpoint.write_text(
        json.dumps(checkpoint_data(merged, classified, "RepeatClassifier"))
    )

    def unexpected_popen(*args, **kwargs):
        raise AssertionError("RepeatClassifier should have been skipped")

    monkeypatch.setattr(reclassify_module.subprocess, "Popen", unexpected_popen)

    result = run_repeatclassifier(
        merged,
        tmp_path,
        "RepeatClassifier",
        logging.getLogger("test"),
    )

    assert classified_output_has_complete_family_set(merged, classified)
    assert result == classified


def test_exact_output_without_completion_checkpoint_is_rerun(
    tmp_path,
    monkeypatch,
):
    merged, classified = write_matching_fastas(tmp_path)

    class FakeProcess:
        def __init__(self):
            self.stdout = []

        def wait(self):
            classified.write_text(">family-1#DNA/TcMar\nAAAA\n")
            return 0

    monkeypatch.setattr(
        reclassify_module.subprocess,
        "Popen",
        lambda *args, **kwargs: FakeProcess(),
    )

    result = run_repeatclassifier(
        merged,
        tmp_path,
        "RepeatClassifier",
        logging.getLogger("test"),
    )

    assert result == classified
    assert (tmp_path / "merged.fa.classified.stale").exists()
    assert (tmp_path / "repeatclassifier.checkpoint.json").exists()


def test_successful_process_is_rejected_when_a_family_length_changes(
    tmp_path,
    monkeypatch,
):
    merged = tmp_path / "merged.fa"
    merged.write_text(">family-1#Unknown\nACGT\n")
    classified = tmp_path / "merged.fa.classified"

    class FakeProcess:
        def __init__(self):
            self.stdout = []

        def wait(self):
            classified.write_text(">family-1#DNA/TcMar\nACG\n")
            return 0

    monkeypatch.setattr(
        reclassify_module.subprocess,
        "Popen",
        lambda *args, **kwargs: FakeProcess(),
    )

    try:
        run_repeatclassifier(
            merged,
            tmp_path,
            "RepeatClassifier",
            logging.getLogger("test"),
        )
    except RuntimeError as error:
        assert "1 changed lengths" in str(error)
    else:
        raise AssertionError("A length-changing output should have been rejected")

    assert json.loads((tmp_path / "repeatclassifier.checkpoint.json").read_text())["status"] == "running"


def test_completed_legacy_manifest_is_migrated(tmp_path, monkeypatch):
    merged, classified = write_matching_fastas(tmp_path)
    manifest = tmp_path / "reclassify.manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "stage": "reclassify",
                "merged_library": str(merged),
                "classified_library": str(classified),
            }
        )
    )

    def unexpected_popen(*args, **kwargs):
        raise AssertionError("Completed legacy output should have been reused")

    monkeypatch.setattr(reclassify_module.subprocess, "Popen", unexpected_popen)

    result = run_repeatclassifier(
        merged,
        tmp_path,
        "RepeatClassifier",
        logging.getLogger("test"),
    )

    assert result == classified
    assert (tmp_path / "repeatclassifier.checkpoint.json").exists()
