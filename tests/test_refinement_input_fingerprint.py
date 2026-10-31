import json
import os
from pathlib import Path

import pandas as pd

from drayte.pipeline.runner import refinement_checkpoint_is_current
from drayte.refinement.refine_repeatmasker import refine_repeatmasker
from drayte.reporting.SummaryFilesGen import required_summary_files, summary_files_are_current


def write_hits(path, repeat_class):
    path.write_text(
        "   SW  perc div. perc del. perc ins. query sequence\n"
        "score  div. del. ins. sequence begin end (left) strand repeat class/family\n"
        f"200 3.0 0.0 0.0 scaffold_1 100 150 (900) + TE_1 {repeat_class} 1 51 (0) 1\n"
    )


def current(result, rmout):
    return refinement_checkpoint_is_current(
        required_outputs=[Path(result[key]) for key in ("refined_tsv", "filtered_tsv", "stats_tsv")],
        manifest_path=Path(result["manifest_json"]),
        expected_parameters={"mode": "loose", "max_gap": 150},
        rmout=rmout,
    )


def test_changed_repeatmasker_class_regenerates_refinement_and_invalidates_summary(tmp_path):
    rmout = tmp_path / "genome.fa.out"
    write_hits(rmout, "Unknown")
    original_mtime = rmout.stat().st_mtime_ns
    result = refine_repeatmasker(rmout, "Test", tmp_path / "refinement", mode="loose", resolve_overlaps=True)
    assert current(result, rmout)
    refined = Path(result["refined_tsv"])
    before = pd.read_csv(refined, sep="\t")
    assert set(before["class"]) == {"Unknown"}

    # Model an existing, current summary before the assisted rerun.
    old = 1_700_000_000_000_000_000
    os.utime(refined, ns=(old, old))
    summaries = tmp_path / "summary"
    summaries.mkdir()
    for path in required_summary_files(summaries, "Test"):
        path.write_text("old summary\n")
        os.utime(path, ns=(old + 1_000_000_000, old + 1_000_000_000))
    assert summary_files_are_current(summaries, "Test", refined)

    # Identity is based on contents, including labels, even with unchanged mtime.
    write_hits(rmout, "RC/Helitron")
    os.utime(rmout, ns=(original_mtime, original_mtime))
    assert not current(result, rmout)
    updated = refine_repeatmasker(rmout, "Test", tmp_path / "refinement", mode="loose", resolve_overlaps=True)
    assert current(updated, rmout)
    assert updated["repeatmasker_input_sha256"] != result["repeatmasker_input_sha256"]
    assert set(pd.read_csv(updated["refined_tsv"], sep="\t")["class"]) == {"RC/Helitron"}
    assert not summary_files_are_current(summaries, "Test", Path(updated["refined_tsv"]))


def test_refinement_reuses_identical_input_despite_changed_timestamp(tmp_path):
    rmout = tmp_path / "genome.fa.out"
    write_hits(rmout, "RC/Helitron")
    result = refine_repeatmasker(rmout, "Test", tmp_path / "refinement", mode="loose", resolve_overlaps=True)
    os.utime(rmout, ns=(1_700_000_000_000_000_000, 1_700_000_000_000_000_000))
    assert current(result, rmout)


def test_legacy_refinement_without_input_hash_is_rebuilt_once(tmp_path):
    rmout = tmp_path / "genome.fa.out"
    write_hits(rmout, "Unknown")
    result = refine_repeatmasker(rmout, "Test", tmp_path / "refinement", mode="loose", resolve_overlaps=True)
    manifest = Path(result["manifest_json"])
    legacy = json.loads(manifest.read_text())
    legacy.pop("repeatmasker_input_sha256")
    manifest.write_text(json.dumps(legacy))
    assert not current(result, rmout)
