import json
import logging
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from drayte import step1_repmodannotation as step1
from drayte.pipeline import final_annotation
from drayte.storage.archive import ArchiveStore
from drayte.utils.repeatmasker_checkpoint import checkpoint_path, legacy_outputs_usable


@pytest.fixture(params=["initial", "final"])
def stage(tmp_path, request):
    root = tmp_path / "run"
    genome = root / "discovery/assemblies_dir/T.fa"
    genome.parent.mkdir(parents=True)
    genome.write_text(">chr1\nACGTACGT\n")
    library = tmp_path / "library.fa"
    library.write_text(">family#RC/Helitron\nACGT\n")
    calls = tmp_path / "calls"
    fail = tmp_path / "fail-once"
    omit = tmp_path / "omit-out"
    script = tmp_path / "RepeatMasker"
    script.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        f"calls = Path({str(calls)!r})\n"
        "with calls.open('a') as out: out.write('called\\n')\n"
        "genome = Path(sys.argv[-1])\n"
        "for suffix in ['.out', '.tbl', '.masked', '.cat.gz', '.align', '.gff']:\n"
        f"    if suffix == '.out' and Path({str(omit)!r}).exists(): continue\n"
        "    Path(str(genome) + suffix).write_text('fresh ' + suffix + '\\n')\n"
        f"fail = Path({str(fail)!r})\n"
        "if fail.exists():\n"
        "    fail.unlink()\n"
        "    sys.exit(122)\n"
    )
    script.chmod(0o755)
    config = SimpleNamespace(outdir_path=root, species="T", threads=2,
                             extra={"repeatmasker_bin": str(script)})
    logger = logging.getLogger("test.repeatmasker.checkpoint")
    if request.param == "initial":
        outdir = root / "discovery/rmasker_dir"

        def run():
            return step1.run_step1(genome=genome, outdir=root / "discovery", species="T",
                                  threads=2, repeatmodeler_dir=tmp_path / "unused-rm",
                                  repeatscout_dir=tmp_path / "unused-rs",
                                  repeatmasker_bin=str(script), repeatmodeler_library=library,
                                  logger=logger)
    else:
        outdir = root / "final_annotation"

        def run():
            return final_annotation.run(config, {"classified_library": str(library)}, logger)

    return SimpleNamespace(run=run, root=root, outdir=outdir, library=library, genome=genome,
                           calls=calls, fail=fail, omit=omit, name=request.param)


def number_calls(stage):
    return stage.calls.read_text().count("called") if stage.calls.exists() else 0


def test_successful_checkpoint_ignores_unchanged_input_timestamps(stage):
    stage.run()
    checkpoint = checkpoint_path(stage.outdir)
    record = json.loads(checkpoint.read_text())
    assert record["state"] == "complete"
    assert not record["legacy_adopted"]
    original_checkpoint = checkpoint.read_bytes()
    stage.library.write_text(stage.library.read_text())
    stage.genome.write_text(stage.genome.read_text())
    stage.run()
    assert number_calls(stage) == 1
    assert checkpoint.read_bytes() == original_checkpoint


def test_changed_library_content_reruns_repeatmasker(stage):
    stage.run()
    stage.library.write_text(">family#LINE/L1\nAAAACCCC\n")
    stage.run()
    assert number_calls(stage) == 2


def test_archived_repeatmasker_products_restore_without_rerunning(stage):
    stage.run()
    checkpoint = checkpoint_path(stage.outdir)
    original_checkpoint = checkpoint.read_bytes()
    products = [stage.outdir / "T.fa.masked"]
    if stage.name == "final":
        products.append(stage.outdir / "T.fa.align")
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in products}
    store = ArchiveStore(stage.root)
    for path in products:
        store.archive(path.relative_to(stage.root))
        assert not path.exists()

    store.restore_all()
    for path, expected in before.items():
        assert (path.read_bytes(), path.stat().st_mtime_ns) == expected
    # Resuming curation/classification can rewrite an unchanged library.
    # Content, rather than its newer timestamp, must control RepeatMasker reuse.
    stage.library.write_text(stage.library.read_text())
    stage.run()
    assert number_calls(stage) == 1
    assert checkpoint.read_bytes() == original_checkpoint


def test_changed_required_output_reruns_repeatmasker(stage):
    stage.run()
    (stage.outdir / "T.fa.out").write_text("tampered output\n")
    stage.run()
    assert number_calls(stage) == 2
    assert (stage.outdir / "T.fa.out").read_text() == "fresh .out\n"


def test_failed_attempt_even_with_complete_products_cannot_be_adopted(stage):
    stage.fail.touch()
    with pytest.raises(RuntimeError, match="exit code 122"):
        stage.run()
    assert json.loads(checkpoint_path(stage.outdir).read_text())["state"] == "running"
    assert Path(str(stage.genome) + ".out").exists()
    stage.run()
    assert number_calls(stage) == 2
    assert json.loads(checkpoint_path(stage.outdir).read_text())["state"] == "complete"


def test_table_only_partial_run_is_not_reused(stage):
    Path(str(stage.genome) + ".tbl").write_text("partial table\n")
    stage.outdir.mkdir(parents=True, exist_ok=True)
    (stage.outdir / "T.fa.tbl").write_text("partial published table\n")
    stage.run()
    assert number_calls(stage) == 1
    assert (stage.outdir / "T.fa.tbl").read_text() == "fresh .tbl\n"


def test_zero_exit_without_required_output_is_not_completion(stage):
    stage.omit.touch()
    # Old published outputs must not satisfy validation of a new attempt.
    stage.outdir.mkdir(parents=True, exist_ok=True)
    (stage.outdir / "T.fa.out").write_text("old output\n")
    with pytest.raises(RuntimeError, match="did not produce required"):
        stage.run()
    assert json.loads(checkpoint_path(stage.outdir).read_text())["state"] == "running"
    stage.omit.unlink()
    stage.run()
    assert number_calls(stage) == 2


def test_complete_legacy_outputs_with_latest_failed_log_are_not_adopted(stage):
    stage.outdir.mkdir(parents=True, exist_ok=True)
    for suffix in (".out", ".tbl", ".masked", ".cat.gz", ".align"):
        (stage.outdir / f"T.fa{suffix}").write_text("legacy output\n")
    (stage.outdir / "T.RepeatMasker.log").write_text(
        "=== RepeatMasker started old ===\n=== RepeatMasker completed successfully ===\n"
        "=== RepeatMasker started new ===\n=== RepeatMasker failed (exit code 122); see full log ===\n"
    )
    stage.run()
    assert number_calls(stage) == 1
    assert not json.loads(checkpoint_path(stage.outdir).read_text())["legacy_adopted"]


def test_legacy_log_started_but_not_completed_is_not_success(tmp_path):
    output = tmp_path / "T.fa.out"
    output.write_text("nonempty\n")
    log = tmp_path / "T.RepeatMasker.log"
    log.write_text("=== RepeatMasker started killed-job ===\n")
    assert not legacy_outputs_usable([output], log)


def test_initial_legacy_source_pair_is_moved_before_checkpoint(tmp_path, stage):
    if stage.name != "initial":
        pytest.skip("discovery-only source-output migration")
    for suffix in (".out", ".tbl"):
        Path(str(stage.genome) + suffix).write_text(f"legacy {suffix}\n")
    stage.run()
    assert number_calls(stage) == 0
    assert not Path(str(stage.genome) + ".out").exists()
    assert (stage.outdir / "T.fa.out").read_text() == "legacy .out\n"
    assert json.loads(checkpoint_path(stage.outdir).read_text())["legacy_adopted"]
