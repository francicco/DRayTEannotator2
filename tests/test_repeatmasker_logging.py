import json
import logging
import os
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

import pytest

from drayte import step1_repmodannotation as step1
from drayte.pipeline import discovery, final_annotation
from drayte.utils import subprocess as command_utils


STDOUT_PAYLOAD = "raw-repeatmasker-stdout-payload"
STDERR_PAYLOAD = "raw-repeatmasker-stderr-payload"
SUFFIXES = [".out", ".tbl", ".gff", ".masked", ".cat.gz", ".align"]


def make_script(tmp_path, name="fake-repeatmasker", *, fail=False, outputs=False):
    script = tmp_path / name
    script.write_text(
        f"#!{sys.executable}\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        f"print({STDOUT_PAYLOAD!r})\n"
        f"print({STDERR_PAYLOAD!r}, file=sys.stderr)\n"
        "print('ARGS=' + json.dumps(sys.argv[1:]))\n"
        "print('CWD=' + os.getcwd())\n"
        "print('MARKER=' + os.environ.get('DRAYTE_LOG_TEST', 'unset'))\n"
        + ("sys.exit(7)\n" if fail else "")
        + (
            "genome = Path(sys.argv[-1])\n"
            f"for suffix in {SUFFIXES!r}:\n"
            "    Path(str(genome) + suffix).write_text('output ' + suffix + '\\n')\n"
            if outputs
            else ""
        )
    )
    script.chmod(0o755)
    return script


def assert_raw_output_is_quiet(caplog, capfd):
    captured = capfd.readouterr()
    for payload in [STDOUT_PAYLOAD, STDERR_PAYLOAD]:
        assert payload not in captured.out
        assert payload not in captured.err
        assert payload not in caplog.text


def test_file_logging_captures_both_streams_and_appends_attempts(tmp_path, caplog, capfd):
    script = make_script(tmp_path)
    logfile = tmp_path / "nested logs" / "T.RepeatMasker.log"
    logfile.parent.mkdir()
    logfile.write_text("previous attempt\n")
    logger = logging.getLogger("test.repeatmasker.helper")
    caplog.set_level(logging.INFO)
    env = {**os.environ, "DRAYTE_LOG_TEST": "environment-was-passed"}

    for _ in range(2):
        command_utils.run_command_to_file(
            [sys.executable, str(script), "argument with spaces"],
            logfile,
            logger,
            cwd=tmp_path,
            env=env,
            label="RepeatMasker",
        )

    text = logfile.read_text()
    assert text.startswith("previous attempt\n")
    assert text.count(STDOUT_PAYLOAD) == 2
    assert text.count(STDERR_PAYLOAD) == 2
    assert text.count("environment-was-passed") == 2
    assert text.count(f"CWD={tmp_path}") == 2
    assert text.count(str(script)) == 2
    assert len(re.findall(r"\d{4}-\d{2}-\d{2}", text)) >= 2
    assert "RepeatMasker" in caplog.text
    assert str(logfile) in caplog.text
    assert_raw_output_is_quiet(caplog, capfd)


def test_failed_command_keeps_diagnostics_in_file_and_reports_log(tmp_path, caplog, capfd):
    script = make_script(tmp_path, fail=True)
    logfile = tmp_path / "T.RepeatMasker.log"
    caplog.set_level(logging.INFO)

    with pytest.raises(RuntimeError) as raised:
        command_utils.run_command_to_file(
            [sys.executable, str(script)],
            logfile,
            logging.getLogger("test.repeatmasker.failure"),
            label="RepeatMasker",
        )

    assert str(logfile) in str(raised.value)
    assert isinstance(raised.value.__cause__, subprocess.CalledProcessError)
    assert raised.value.__cause__.returncode == 7
    assert STDOUT_PAYLOAD in logfile.read_text()
    assert STDERR_PAYLOAD in logfile.read_text()
    errors = [record.getMessage() for record in caplog.records if record.levelno >= logging.ERROR]
    assert any(str(logfile) in error for error in errors)
    assert_raw_output_is_quiet(caplog, capfd)


def test_missing_executable_reports_attempt_log_and_chained_error(tmp_path, caplog, capfd):
    executable = tmp_path / "missing-repeatmasker"
    logfile = tmp_path / "T.RepeatMasker.log"
    caplog.set_level(logging.INFO)

    with pytest.raises(RuntimeError) as raised:
        command_utils.run_command_to_file(
            [str(executable)],
            logfile,
            logging.getLogger("test.repeatmasker.missing"),
            label="RepeatMasker",
        )

    assert isinstance(raised.value.__cause__, FileNotFoundError)
    assert str(logfile) in str(raised.value)
    assert logfile.is_file()
    assert str(executable) in logfile.read_text()
    assert any(record.levelno >= logging.ERROR for record in caplog.records)
    assert_raw_output_is_quiet(caplog, capfd)


@pytest.fixture
def stage_inputs(tmp_path):
    genome = tmp_path / "source.fa"
    genome.write_text(">chr1\nACGTACGT\n")
    library = tmp_path / "library.fa"
    library.write_text(">rnd-1_family-1#LINE/L1\nACGT\n")
    root = tmp_path / "run"
    assemblies = root / "discovery" / "assemblies_dir"
    assemblies.mkdir(parents=True)
    (assemblies / "T.fa").write_text(genome.read_text())
    script = make_script(tmp_path, outputs=True)
    logger = logging.getLogger("test.repeatmasker.stages")
    config = SimpleNamespace(
        outdir_path=root,
        species="T",
        threads=3,
        extra={"repeatmasker_bin": str(script)},
    )

    def initial():
        return step1.run_step1(
            genome=genome,
            outdir=root / "discovery",
            species="T",
            threads=3,
            repeatmodeler_dir=tmp_path / "unused-repeatmodeler",
            repeatscout_dir=tmp_path / "unused-repeatscout",
            repeatmasker_bin=str(config.extra["repeatmasker_bin"]),
            repeatmodeler_library=library,
            logger=logger,
        )

    def final():
        return final_annotation.run(config, {"classified_library": str(library)}, logger)

    return SimpleNamespace(
        root=root,
        genome=assemblies / "T.fa",
        library=library,
        config=config,
        script=script,
        initial=initial,
        final=final,
    )


@pytest.mark.parametrize("stage", ["initial", "final"])
def test_repeatmasker_stage_logs_raw_output_and_preserves_arguments_and_outputs(
    stage_inputs, stage, caplog, capfd
):
    caplog.set_level(logging.INFO)
    result = getattr(stage_inputs, stage)()
    outdir = (
        stage_inputs.root / "discovery" / "rmasker_dir"
        if stage == "initial"
        else stage_inputs.root / "final_annotation"
    )
    logfile = outdir / "T.RepeatMasker.log"
    assert result["repeatmasker_log"] == str(logfile)
    text = logfile.read_text()
    assert STDOUT_PAYLOAD in text
    assert STDERR_PAYLOAD in text
    args = json.loads(next(line[5:] for line in text.splitlines() if line.startswith("ARGS=")))
    assert args[-1] == str(stage_inputs.genome)
    assert args[args.index("-pa") + 1] == "3"
    expected_library = (
        stage_inputs.root / "discovery" / "rmodeler_dir" / "consensi.fa.classified"
        if stage == "initial"
        else stage_inputs.library
    )
    assert args[args.index("-lib") + 1] == str(expected_library)
    assert all(option in args for option in ["-s", "-xsmall", "-gc", "-gff", "-norna"])
    if stage == "final":
        assert "-a" in args
        assert f"CWD={stage_inputs.genome.parent}" in text
    for suffix in SUFFIXES if stage == "final" else SUFFIXES[:-1]:
        assert (outdir / f"T.fa{suffix}").read_text() == f"output {suffix}\n"
        assert not Path(str(stage_inputs.genome) + suffix).exists()
    assert_raw_output_is_quiet(caplog, capfd)

    before = logfile.read_bytes()
    before_mtime = logfile.stat().st_mtime_ns
    resumed = getattr(stage_inputs, stage)()
    assert resumed["repeatmasker_log"] == str(logfile)
    assert logfile.read_bytes() == before
    assert logfile.stat().st_mtime_ns == before_mtime


@pytest.mark.parametrize("stage", ["initial", "final"])
@pytest.mark.parametrize("has_log", [False, True])
def test_cached_stage_does_not_create_or_truncate_repeatmasker_log(stage_inputs, stage, has_log):
    outdir = (
        stage_inputs.root / "discovery" / "rmasker_dir"
        if stage == "initial"
        else stage_inputs.root / "final_annotation"
    )
    outdir.mkdir(parents=True)
    for suffix in SUFFIXES:
        (outdir / f"T.fa{suffix}").write_text(f"cached {suffix}\n")
    logfile = outdir / "T.RepeatMasker.log"
    if has_log:
        logfile.write_text("historical diagnostic log\n")
        before = logfile.stat().st_mtime_ns
    stage_inputs.config.extra["repeatmasker_bin"] = str(outdir / "does-not-exist")

    result = getattr(stage_inputs, stage)()

    if has_log:
        assert result["repeatmasker_log"] == str(logfile)
        assert logfile.read_text() == "historical diagnostic log\n"
        assert logfile.stat().st_mtime_ns == before
    else:
        assert "repeatmasker_log" not in result
        assert not logfile.exists()


@pytest.mark.parametrize("stage", ["initial", "final"])
def test_stage_failure_reports_log_without_echoing_raw_diagnostics(
    tmp_path, stage_inputs, stage, caplog, capfd
):
    stage_inputs.config.extra["repeatmasker_bin"] = str(
        make_script(tmp_path, name="failed-repeatmasker", fail=True)
    )
    caplog.set_level(logging.INFO)
    with pytest.raises(RuntimeError) as raised:
        getattr(stage_inputs, stage)()
    outdir = (
        stage_inputs.root / "discovery" / "rmasker_dir"
        if stage == "initial"
        else stage_inputs.root / "final_annotation"
    )
    logfile = outdir / "T.RepeatMasker.log"
    assert str(logfile) in str(raised.value)
    assert STDOUT_PAYLOAD in logfile.read_text()
    assert STDERR_PAYLOAD in logfile.read_text()
    assert_raw_output_is_quiet(caplog, capfd)


def test_discovery_forwards_initial_repeatmasker_log_to_stage_result(
    tmp_path, stage_inputs, caplog, capfd
):
    config = stage_inputs.config
    config.genome_path = stage_inputs.genome
    config.extra.update(
        repeatmodeler_dir=str(tmp_path / "unused-repeatmodeler"),
        repeatscout_dir=str(tmp_path / "unused-repeatscout"),
        repeatmodeler_library=str(stage_inputs.library),
    )
    caplog.set_level(logging.INFO)

    result = discovery.run(config, logging.getLogger("test.repeatmasker.discovery"))

    logfile = stage_inputs.root / "discovery" / "rmasker_dir" / "T.RepeatMasker.log"
    assert result["stage"] == "discovery"
    assert result["repeatmasker_log"] == str(logfile)
    assert STDOUT_PAYLOAD in logfile.read_text()
    assert STDERR_PAYLOAD in logfile.read_text()
    assert json.loads(json.dumps(result))["repeatmasker_log"] == str(logfile)
    assert_raw_output_is_quiet(caplog, capfd)
