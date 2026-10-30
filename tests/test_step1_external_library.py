import logging
from pathlib import Path

import pytest

from drayte import step1_repmodannotation as step1


def test_external_library_skips_repeatmodeler_and_reuses_common_outputs(
    tmp_path, monkeypatch
):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGTACGT\n")

    external_library = tmp_path / "consensi.fa.classified"
    external_library.write_text(">rnd-1_family-1#LINE/L1\nACGT\n")

    outdir = tmp_path / "discovery"
    stale_dir = outdir / "rmodeler_dir"
    stale_dir.mkdir(parents=True)
    (stale_dir / "consensi.fa.classified").write_text(">stale#Unknown\nAAAA\n")
    (stale_dir / "Test-families.fa").write_text(">stale#Unknown\nAAAA\n")
    (stale_dir / "Test-families.mod.fa").write_text(">stale#Unknown\nAAAA\n")

    commands = []

    def record_command(cmd, cwd=None, logger=None):
        commands.append([str(value) for value in cmd])

    def record_file_command(cmd, log_file, logger, cwd=None, env=None, label="Command"):
        commands.append([str(value) for value in cmd])
        Path(log_file).write_text("recorded RepeatMasker\n")
        for suffix in (".out", ".tbl"):
            Path(str(cmd[-1]) + suffix).write_text(f"RepeatMasker {suffix}\n")

    monkeypatch.setattr(step1, "run_command", record_command)
    monkeypatch.setattr(step1, "run_command_to_file", record_file_command)

    result = step1.run_step1(
        genome=genome,
        outdir=outdir,
        species="Test",
        threads=2,
        repeatmodeler_dir=tmp_path / "unused-repeatmodeler",
        repeatscout_dir=tmp_path / "unused-repeatscout",
        repeatmasker_bin="RepeatMasker",
        repeatmodeler_library=external_library,
        logger=logging.getLogger("test.external-library"),
    )

    assert [Path(command[0]).name for command in commands] == ["RepeatMasker"]
    assert Path(result["raw_library"]).read_text() == external_library.read_text()
    assert Path(result["species_families"]).read_text() == external_library.read_text()
    assert (
        Path(result["normalized_library"]).read_text()
        == ">Test-rnd-1_family-1#LINE/L1\nACGT\n"
    )
    assert Path(result["rmasker_dir"]) == outdir.resolve() / "rmasker_dir"


@pytest.mark.parametrize("empty_cached_library", [False, True])
def test_normal_path_still_runs_builddatabase_and_repeatmodeler(
    tmp_path, monkeypatch, empty_cached_library
):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGTACGT\n")

    repeatmodeler_dir = tmp_path / "RepeatModeler"
    repeatmodeler_dir.mkdir()
    for name in ("BuildDatabase", "RepeatModeler"):
        (repeatmodeler_dir / name).touch()

    repeatscout_dir = tmp_path / "RepeatScout"
    repeatscout_dir.mkdir()
    for name in ("RepeatScout", "build_lmer_table"):
        (repeatscout_dir / name).touch()

    if empty_cached_library:
        cached = tmp_path / "discovery/rmodeler_dir/consensi.fa.classified"
        cached.parent.mkdir(parents=True)
        cached.touch()

    commands = []
    repeatmodeler_commands = []

    def record_command(cmd, cwd=None, logger=None):
        commands.append([str(value) for value in cmd])

    def record_file_command(cmd, log_file, logger, cwd=None, env=None, label="Command"):
        commands.append([str(value) for value in cmd])
        Path(log_file).write_text("recorded RepeatMasker\n")
        for suffix in (".out", ".tbl"):
            Path(str(cmd[-1]) + suffix).write_text(f"RepeatMasker {suffix}\n")

    def fake_repeatmodeler(cmd, log_file, cwd=None, logger=None, prefix="subprocess"):
        repeatmodeler_commands.append([str(value) for value in cmd])
        result_dir = cwd / "RM_test"
        result_dir.mkdir()
        (result_dir / "consensi.fa.classified").write_text(
            ">rnd-2_family-3#DNA/TcMar\nACGT\n"
        )

    monkeypatch.setattr(step1, "run_command", record_command)
    monkeypatch.setattr(step1, "run_command_to_file", record_file_command)
    monkeypatch.setattr(
        step1, "run_command_to_logger_and_file", fake_repeatmodeler
    )

    result = step1.run_step1(
        genome=genome,
        outdir=tmp_path / "discovery",
        species="Test",
        threads=2,
        repeatmodeler_dir=repeatmodeler_dir,
        repeatscout_dir=repeatscout_dir,
        repeatmasker_bin="RepeatMasker",
        logger=logging.getLogger("test.normal-path"),
    )

    assert [Path(command[0]).name for command in commands] == [
        "BuildDatabase",
        "RepeatMasker",
    ]
    assert [Path(command[0]).name for command in repeatmodeler_commands] == [
        "RepeatModeler"
    ]
    assert Path(result["raw_library"]).read_text().startswith(
        ">rnd-2_family-3#DNA/TcMar"
    )


@pytest.mark.parametrize("annotation_cached", [False, True])
def test_native_cached_library_needs_no_database_and_preserves_annotation_resume(
    tmp_path, monkeypatch, annotation_cached
):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGTACGT\n")
    repeatmodeler_dir = tmp_path / "RepeatModeler"
    repeatmodeler_dir.mkdir()
    for name in ("BuildDatabase", "RepeatModeler"):
        (repeatmodeler_dir / name).touch()
    repeatscout_dir = tmp_path / "RepeatScout"
    repeatscout_dir.mkdir()
    for name in ("RepeatScout", "build_lmer_table"):
        (repeatscout_dir / name).touch()

    outdir = tmp_path / "discovery"
    rmodeler_dir = outdir / "rmodeler_dir"
    rmodeler_dir.mkdir(parents=True)
    library = rmodeler_dir / "consensi.fa.classified"
    library.write_text(">rnd-1_family-1#LINE/L1\nACGT\n")
    rmasker_dir = outdir / "rmasker_dir"
    rmasker_dir.mkdir()
    cached_bytes = {}
    cached_mtimes = {}
    if annotation_cached:
        for suffix in (".out", ".tbl"):
            path = rmasker_dir / f"Test.fa{suffix}"
            path.write_text(f"cached RepeatMasker {suffix}\n")
            cached_bytes[suffix] = path.read_bytes()
            cached_mtimes[suffix] = path.stat().st_mtime_ns

    def unexpected_discovery(*args, **kwargs):
        raise AssertionError("Cached native library must skip BuildDatabase and RepeatModeler")

    repeatmasker_calls = []

    def fake_repeatmasker(cmd, log_file, logger, **kwargs):
        repeatmasker_calls.append([str(value) for value in cmd])
        Path(log_file).write_text("recorded RepeatMasker\n")
        for suffix in (".out", ".tbl"):
            Path(str(cmd[-1]) + suffix).write_text(f"new RepeatMasker {suffix}\n")

    monkeypatch.setattr(step1, "run_command", unexpected_discovery)
    monkeypatch.setattr(step1, "run_command_to_logger_and_file", unexpected_discovery)
    monkeypatch.setattr(step1, "run_command_to_file", fake_repeatmasker)

    def run():
        return step1.run_step1(
            genome=genome, outdir=outdir, species="Test", threads=2,
            repeatmodeler_dir=repeatmodeler_dir, repeatscout_dir=repeatscout_dir,
            repeatmasker_bin="RepeatMasker", logger=logging.getLogger("test.cached-native"),
        )

    result = run()
    assert not step1.repeatmodeler_db_exists(outdir / "assemblies_dir", "Test")
    assert Path(result["raw_library"]).read_bytes() == library.read_bytes()
    assert Path(result["normalized_library"]).read_text() == ">Test-rnd-1_family-1#LINE/L1\nACGT\n"
    assert len(repeatmasker_calls) == (0 if annotation_cached else 1)
    if annotation_cached:
        for suffix in (".out", ".tbl"):
            path = rmasker_dir / f"Test.fa{suffix}"
            assert path.read_bytes() == cached_bytes[suffix]
            assert path.stat().st_mtime_ns == cached_mtimes[suffix]
    # A second resume must still skip discovery and reuse the validated
    # RepeatMasker completion record, despite the database remaining absent.
    calls_before = len(repeatmasker_calls)
    checkpoint = rmasker_dir / "repeatmasker.checkpoint.json"
    checkpoint_bytes = checkpoint.read_bytes()
    run()
    assert len(repeatmasker_calls) == calls_before
    assert checkpoint.read_bytes() == checkpoint_bytes
