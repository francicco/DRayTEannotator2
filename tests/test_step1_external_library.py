import logging
from pathlib import Path

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

    monkeypatch.setattr(step1, "run_command", record_command)

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


def test_normal_path_still_runs_builddatabase_and_repeatmodeler(tmp_path, monkeypatch):
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

    commands = []
    repeatmodeler_commands = []

    def record_command(cmd, cwd=None, logger=None):
        commands.append([str(value) for value in cmd])

    def fake_repeatmodeler(cmd, log_file, cwd=None, logger=None, prefix="subprocess"):
        repeatmodeler_commands.append([str(value) for value in cmd])
        result_dir = cwd / "RM_test"
        result_dir.mkdir()
        (result_dir / "consensi.fa.classified").write_text(
            ">rnd-2_family-3#DNA/TcMar\nACGT\n"
        )

    monkeypatch.setattr(step1, "run_command", record_command)
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
