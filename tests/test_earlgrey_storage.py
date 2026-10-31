import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from drayte.storage.archive import ArchiveStore
from drayte.storage.policy import output_lock


SCRIPT = Path(__file__).resolve().parents[1] / "compact_earlgrey.py"
spec = importlib.util.spec_from_file_location("standalone_compact_earlgrey", SCRIPT)
earlgrey = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = earlgrey
spec.loader.exec_module(earlgrey)


@pytest.fixture
def run_tree(tmp_path):
    outer = tmp_path / "Test_Earlgrey"
    inner = outer / "Test_EarlGrey"
    inner.mkdir(parents=True)

    def write(relative, text):
        path = inner / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    write("Test_EarlGrey.log", "<<< Done in 01:02:03 >>>\n<<< " + earlgrey.FINAL_MESSAGE + " path >>>\n")
    write("Test_Database/Test-families.fa", ">family\nACGT\n")
    library = write("Test_strainer/Test-families.fa.strained", ">family\nACGTACGT\n")
    write("Test_summaryFiles/" + library.name, library.read_text())
    for suffix in ("highLevelCount.txt", "highLevelCount.kable", "summaryPie.pdf"):
        write("Test_summaryFiles/Test." + suffix, "complete\n")
    write("Test_summaryFiles/Test_divergence_summary_table.tsv", "complete\n")
    for suffix in ("bed", "gff"):
        write("Test_mergedRepeats/looseMerge/Test.filteredRepeats." + suffix, "annotation\n")
        write("Test_summaryFiles/Test.filteredRepeats." + suffix, "annotation\n")
    for suffix in ("out", "tbl"):
        write("Test_RepeatMasker_Against_Custom_Library/genome.fa.prep." + suffix, "hits\n")
    for suffix in ("align", "masked", "cat.gz"):
        write("Test_RepeatMasker_Against_Custom_Library/genome.fa.prep." + suffix, "AAAA" * 10000)
    write("Test_RepeatLandscape/landscape.dat", "AAAA" * 10000)
    write("Test_heliano/HEL_1/RC.representative.gff", "annotation\n")
    write("Test_heliano/HEL_1/HLE1/scratch.fa", "AAAA" * 10000)
    for index in range(25):
        write(f"Test_RepeatModeler/RM_old/chunk_{index}.fa", "ACGT" * 1000)
        write(f"Test_strainer/TS_old/chunk_{index}.fa", "ACGT" * 1000)
    return outer, inner, write


def tree_files(root):
    return {str(path.relative_to(root)): path.read_bytes()
            for path in root.rglob("*") if path.is_file()}


def test_preview_is_read_only_and_accepts_outer_or_inner(run_tree):
    outer, inner, _ = run_tree
    before = tree_files(outer)
    assert earlgrey.layout(outer) == (outer, inner, "Test")
    assert earlgrey.layout(inner) == (outer, inner, "Test")
    report = earlgrey.preview(inner, "Test")
    assert report["eligible"]
    assert report["read_only"]
    assert any(item["path"] == "Test_RepeatLandscape" for item in report["plan"])
    assert tree_files(outer) == before
    assert not (inner / ".storage").exists()


def test_single_copied_script_compacts_and_restores_without_installed_packages(run_tree, tmp_path):
    outer, inner, _ = run_tree
    isolated = tmp_path / "independent_script"
    isolated.mkdir()
    copied = isolated / SCRIPT.name
    shutil.copy2(SCRIPT, copied)
    assert list(isolated.iterdir()) == [copied]
    before = tree_files(inner)
    mtimes = {relative: (inner / relative).stat().st_mtime_ns for relative in before}

    def run(*options):
        # -I removes the repository/PYTHONPATH; -S disables site packages. Only
        # the single copied file and Python's standard library are available.
        result = subprocess.run(
            [sys.executable, "-I", "-S", str(copied), "--root", str(outer), *options],
            cwd=isolated, check=True, capture_output=True, text=True,
        )
        return json.loads(result.stdout)

    preview = run()
    assert preview["eligible"] and preview["read_only"]
    assert tree_files(inner) == before
    assert not (inner / ".storage").exists()
    report = run("--apply", "--idle")
    assert report["actions"]
    assert report["after"]["allocated_bytes"] < report["before"]["allocated_bytes"]
    for relative, data in before.items():
        if "summaryFiles/" in relative or "Database/" in relative or relative.endswith((".out", ".tbl", ".cat.gz", ".strained")):
            assert (inner / relative).read_bytes() == data
    catalog = json.loads((inner / ".storage/archives.json").read_text())
    for record in catalog["archives"].values():
        assert hashlib.sha256((inner / record["archive"]).read_bytes()).hexdigest() == record["archive_sha256"]
    assert run("--restore")["read_only"]
    restored = run("--restore", "--apply", "--idle")
    assert restored["restored"]
    for relative, data in before.items():
        assert (inner / relative).read_bytes() == data
        assert (inner / relative).stat().st_mtime_ns == mtimes[relative]
    assert not list((inner / ".storage/archives").iterdir())


@pytest.mark.parametrize("creator,restorer", [
    (ArchiveStore, "standalone"), ("standalone", ArchiveStore),
])
def test_existing_catalogs_are_compatible_in_both_directions(run_tree, creator, restorer):
    _, inner, _ = run_tree
    creator = earlgrey.ArchiveStore if creator == "standalone" else creator
    restorer = earlgrey.ArchiveStore if restorer == "standalone" else restorer
    target = "Test_RepeatModeler/RM_old"
    source = inner / target
    before = tree_files(source)
    with output_lock(inner):
        record = creator(inner).archive(target)
        assert not source.exists()
        assert restorer(inner).restore_all() == [target]
    assert tree_files(source) == before
    assert not (inner / record["archive"]).exists()


def test_compact_retains_libraries_final_products_and_restores_exactly(run_tree):
    outer, inner, _ = run_tree
    before = tree_files(inner)
    with output_lock(inner):
        report = earlgrey.compact(inner, "Test")
    assert report["after"]["files_and_symlinks"] < report["before"]["files_and_symlinks"]
    assert report["after"]["allocated_bytes"] < report["before"]["allocated_bytes"]
    for relative in before:
        if "summaryFiles/" in relative or "Database/" in relative or relative.endswith((".out", ".tbl", ".log", ".strained", ".cat.gz")):
            assert (inner / relative).read_bytes() == before[relative]
    assert not (inner / "Test_mergedRepeats").exists()
    assert (inner / "Test_heliano/HEL_1/RC.representative.gff").is_file()
    assert earlgrey.preview(inner, "Test")["eligible"]
    with output_lock(inner):
        assert earlgrey.compact(inner, "Test")["actions"] == []
        ArchiveStore(inner).restore_all()
    for relative, data in before.items():
        assert (inner / relative).read_bytes() == data
    assert earlgrey.validate_completed(inner, "Test")


@pytest.mark.parametrize("relative", ["Test_EarlGrey.log", "Test_summaryFiles/Test.filteredRepeats.bed",
                                      "Test_Database/Test-families.fa",
                                      "Test_RepeatMasker_Against_Custom_Library/genome.fa.prep.tbl"])
def test_incomplete_runs_are_blocked_without_mutating(run_tree, relative):
    outer, inner, _ = run_tree
    (inner / relative).unlink()
    before = tree_files(inner)
    assert not earlgrey.preview(inner, "Test")["eligible"]
    with pytest.raises(ValueError):
        earlgrey.compact(inner, "Test")
    assert tree_files(inner) == before


def test_changed_export_blocks_compaction(run_tree):
    _, inner, write = run_tree
    write("Test_summaryFiles/Test.filteredRepeats.bed", "different\n")
    assert not earlgrey.preview(inner, "Test")["eligible"]
    assert (inner / "Test_mergedRepeats").is_dir()


def test_strict_fallback_and_recovered_early_error_allowed(run_tree):
    _, inner, write = run_tree
    for suffix in ("bed", "gff"):
        source = inner / f"Test_mergedRepeats/looseMerge/Test.filteredRepeats.{suffix}"
        source.rename(source.parent.parent / source.name)
    write("Test_EarlGrey.log", "ERROR: loose merge defragmentation failed, trying strict merge...\n<<< Done in 00:01:00 >>>\n" + earlgrey.FINAL_MESSAGE)
    assert earlgrey.preview(inner, "Test")["eligible"]
    with (inner / "Test_EarlGrey.log").open("a") as handle:
        handle.write("\nERROR: later failure\n")
    assert not earlgrey.preview(inner, "Test")["eligible"]


def test_symlinks_into_archives_block_before_any_changes(run_tree):
    _, inner, _ = run_tree
    (inner / "retained").symlink_to(inner / "Test_RepeatLandscape/landscape.dat")
    before = tree_files(inner)
    assert not earlgrey.preview(inner, "Test")["eligible"]
    with pytest.raises(ValueError, match="retained symlink"):
        earlgrey.compact(inner, "Test")
    assert tree_files(inner) == before


def test_outside_symlink_in_later_candidate_blocks_all_archives(run_tree, tmp_path):
    _, inner, _ = run_tree
    external = tmp_path / "external"
    external.write_text("external")
    (inner / "Test_RepeatLandscape/outside").symlink_to(external)
    with pytest.raises(ValueError, match="escapes archive root"):
        earlgrey.compact(inner, "Test")
    assert (inner / "Test_RepeatModeler/RM_old").is_dir()
    assert not (inner / ".storage").exists()


def test_cli_apply_requires_idle_confirmation_and_lock(run_tree, capsys):
    outer, inner, _ = run_tree
    with pytest.raises(SystemExit) as error:
        earlgrey.main(["--root", str(outer), "--apply"])
    assert error.value.code == 2
    assert "requires --idle" in capsys.readouterr().err
    with output_lock(inner):
        with pytest.raises(SystemExit) as error:
            earlgrey.main(["--root", str(outer), "--apply", "--idle"])
    assert error.value.code == 2
    assert (inner / "Test_RepeatModeler/RM_old").is_dir()


def command(outer):
    return ["earlGrey", "-g", "genome.fa", "-s", "Test", "-o", str(outer), "-t", "2"]


def test_wrapper_restores_before_launch_then_compacts(run_tree, monkeypatch):
    outer, inner, _ = run_tree
    with output_lock(inner):
        earlgrey.compact(inner, "Test")
    calls = []

    def fake_run(args):
        calls.append(args)
        assert (inner / "Test_mergedRepeats/looseMerge/Test.filteredRepeats.bed").is_file()
        assert not ArchiveStore(inner).records()
        return 0

    monkeypatch.setattr(earlgrey, "run_native", fake_run)
    result = earlgrey.run_wrapped(outer, inner, "Test", command(outer))
    assert calls == [command(outer)]
    assert result["compacted"]
    assert not (inner / "Test_mergedRepeats").exists()


def test_wrapper_failure_keeps_working_files(run_tree, monkeypatch):
    outer, inner, _ = run_tree
    monkeypatch.setattr(earlgrey, "run_native", lambda *args, **kwargs: 7)
    result = earlgrey.run_wrapped(outer, inner, "Test", command(outer))
    assert result["returncode"] == 7
    assert not result["compacted"]
    assert (inner / "Test_RepeatModeler/RM_old").is_dir()
    assert not ArchiveStore(inner).records()


def test_wrapper_never_runs_wrong_command_or_wrong_output(run_tree, monkeypatch):
    outer, inner, _ = run_tree
    monkeypatch.setattr(earlgrey, "run_native", lambda *args, **kwargs: pytest.fail("Unexpected tool run"))
    for args in (["other", "-s", "Test", "-o", str(outer)],
                 ["earlGrey", "-g", "genome.fa", "-s", "Other", "-o", str(outer)],
                 ["earlGrey", "-g", "genome.fa", "-s", "Test", "-o", str(outer.parent)]):
        with pytest.raises(ValueError):
            earlgrey.run_wrapped(outer, inner, "Test", args)
    assert not (inner / ".storage").exists()


def test_corrupt_archive_blocks_restore_and_wrapper_launch(run_tree, monkeypatch):
    outer, inner, _ = run_tree
    with output_lock(inner):
        earlgrey.compact(inner, "Test")
    record = next(iter(ArchiveStore(inner).records().values()))
    (inner / record["archive"]).write_bytes(b"corrupt")
    monkeypatch.setattr(earlgrey, "run_native", lambda *args, **kwargs: pytest.fail("Unexpected tool run"))
    with pytest.raises(ValueError, match="checksum"):
        earlgrey.run_wrapped(outer, inner, "Test", command(outer))


def test_cli_restore_preview_and_apply(run_tree, capsys):
    outer, inner, _ = run_tree
    with output_lock(inner):
        earlgrey.compact(inner, "Test")
    assert earlgrey.main(["--root", str(outer), "--restore"]) == 0
    assert json.loads(capsys.readouterr().out)["read_only"]
    assert not (inner / "Test_mergedRepeats").exists()
    assert earlgrey.main(["--root", str(outer), "--restore", "--apply", "--idle"]) == 0
    assert (inner / "Test_mergedRepeats").is_dir()


def test_symlink_or_ambiguous_root_rejected(run_tree, tmp_path):
    outer, inner, _ = run_tree
    link = tmp_path / "Alias_Earlgrey"
    link.symlink_to(outer)
    with pytest.raises(ValueError, match="symlinked"):
        earlgrey.layout(link, "Test")
    (inner / "Test_EarlGrey").mkdir()
    with pytest.raises(ValueError, match="Ambiguous"):
        earlgrey.layout(inner)


def test_empty_precreated_outer_root_is_blocked_read_only(tmp_path, capsys):
    outer = tmp_path / "Test_Earlgrey"
    outer.mkdir()
    assert earlgrey.main(["--root", str(outer)]) == 0
    report = json.loads(capsys.readouterr().out)
    assert not report["eligible"]
    assert not list(outer.iterdir())


def test_wrapper_rejects_combined_duplicate_output_flag(run_tree, monkeypatch):
    outer, inner, _ = run_tree
    monkeypatch.setattr(earlgrey, "run_native", lambda *args, **kwargs: pytest.fail("Unexpected tool run"))
    with pytest.raises(ValueError, match="combined option"):
        earlgrey.run_wrapped(outer, inner, "Test", command(outer) + ["-oelsewhere"])
    assert not (inner / ".storage").exists()


def test_symlink_report_blocks_before_archiving(run_tree, tmp_path):
    _, inner, _ = run_tree
    storage = inner / ".storage"
    storage.mkdir()
    external = tmp_path / "report"
    external.write_text("preserve")
    (storage / "earlgrey_report.json").symlink_to(external)
    with pytest.raises(ValueError, match="symlinked EarlGrey storage report"):
        earlgrey.compact(inner, "Test")
    assert (inner / "Test_RepeatModeler/RM_old").is_dir()
    assert external.read_text() == "preserve"


def test_cli_wrapper_existing_run_requires_idle_before_restore_or_launch(run_tree, monkeypatch, capsys):
    outer, inner, _ = run_tree
    with output_lock(inner):
        earlgrey.compact(inner, "Test")
    monkeypatch.setattr(earlgrey, "run_native", lambda *args, **kwargs: pytest.fail("Unexpected tool run"))
    with pytest.raises(SystemExit) as error:
        earlgrey.main(["--root", str(outer), "--run", *command(outer)])
    assert error.value.code == 2
    assert "requires --idle" in capsys.readouterr().err
    assert not (inner / "Test_mergedRepeats").exists()


def test_preview_catches_corrupt_existing_archive(run_tree):
    _, inner, _ = run_tree
    with output_lock(inner):
        earlgrey.compact(inner, "Test")
    record = next(iter(ArchiveStore(inner).records().values()))
    (inner / record["archive"]).write_bytes(b"corrupt")
    report = earlgrey.preview(inner, "Test")
    assert not report["eligible"]
    assert "checksum" in report["errors"][0]


def test_appended_incomplete_rerun_log_rejects_old_completion(run_tree):
    _, inner, _ = run_tree
    with (inner / "Test_EarlGrey.log").open("a") as handle:
        handle.write("<<< Cleaning Genome >>>\npartial progress\n")
    report = earlgrey.preview(inner, "Test")
    assert not report["eligible"]
    assert "newer stage start" in report["errors"][0]


def test_native_interrupt_stops_owned_group_preserving_exception(monkeypatch):
    import signal
    calls = []

    class Process:
        pid = 765432

        def wait(self, timeout=None):
            calls.append(("wait", timeout))
            if timeout is None:
                raise KeyboardInterrupt("interrupted")
            return -15

    process = Process()

    def start(args, start_new_session):
        assert start_new_session is True
        return process

    monkeypatch.setattr(earlgrey.subprocess, "Popen", start)
    monkeypatch.setattr(earlgrey.os, "killpg", lambda pid, sig: calls.append((pid, sig)))
    with pytest.raises(KeyboardInterrupt, match="interrupted"):
        earlgrey.run_native(["earlGrey"])
    assert (process.pid, signal.SIGTERM) in calls
    assert (process.pid, signal.SIGKILL) in calls
    assert ("wait", 10) in calls


def test_native_failure_escalates_stubborn_group(monkeypatch):
    import signal
    import subprocess
    calls = []

    class Process:
        pid = 765432
        waits = 0

        def wait(self, timeout=None):
            self.waits += 1
            if timeout is None:
                return 7
            if self.waits == 2:
                raise subprocess.TimeoutExpired("earlGrey", timeout)
            return -9

    process = Process()
    monkeypatch.setattr(earlgrey.subprocess, "Popen", lambda *args, **kwargs: process)
    monkeypatch.setattr(earlgrey.os, "killpg", lambda pid, sig: calls.append((pid, sig)))
    assert earlgrey.run_native(["earlGrey"]) == 7
    assert calls == [(process.pid, signal.SIGTERM), (process.pid, signal.SIGKILL)]


def test_native_cleanup_never_signals_parent_group(monkeypatch):
    process = SimpleNamespace(pid=earlgrey.os.getpgrp(), wait=lambda timeout=None: 7)
    monkeypatch.setattr(earlgrey.subprocess, "Popen", lambda *args, **kwargs: process)
    monkeypatch.setattr(earlgrey.os, "killpg", lambda *args: pytest.fail("Signaled parent group"))
    assert earlgrey.run_native(["earlGrey"]) == 7


def test_wrapper_holds_lock_until_interrupt_cleanup(run_tree, monkeypatch):
    outer, inner, _ = run_tree

    def interrupted(args):
        with pytest.raises(RuntimeError, match="holds the output lock"):
            with output_lock(inner):
                pass
        raise KeyboardInterrupt()

    monkeypatch.setattr(earlgrey, "run_native", interrupted)
    with pytest.raises(KeyboardInterrupt):
        earlgrey.run_wrapped(outer, inner, "Test", command(outer))
    with output_lock(inner):
        pass
    assert (inner / "Test_RepeatModeler/RM_old").is_dir()


def test_real_native_failure_stops_descendant(tmp_path):
    import os
    import shutil
    import signal
    import subprocess
    import sys
    import time
    if os.name != "posix" or not shutil.which("ps"):
        pytest.skip("POSIX process inspection is required")
    try:
        subprocess.run(["ps", "-o", "stat=", "-p", str(os.getpid())],
                       capture_output=True, text=True, check=True)
    except (OSError, subprocess.SubprocessError):
        pytest.skip("Process inspection is unavailable in this execution environment")
    marker = tmp_path / "child.pid"
    program = (
        "import subprocess, sys; from pathlib import Path; "
        "child=subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)']); "
        "Path(sys.argv[1]).write_text(str(child.pid)); sys.exit(7)"
    )
    try:
        assert earlgrey.run_native([sys.executable, "-c", program, str(marker)]) == 7
        child_pid = marker.read_text().strip()
        deadline = time.monotonic() + 2
        while True:
            status = subprocess.run(["ps", "-o", "stat=", "-p", child_pid],
                                    capture_output=True, text=True).stdout.strip()
            if not status or status.startswith("Z"):
                break
            if time.monotonic() > deadline:
                pytest.fail(f"Native descendant survived cleanup: {child_pid}: {status}")
            time.sleep(0.02)
    finally:
        if marker.exists():
            try:
                os.kill(int(marker.read_text()), signal.SIGKILL)
            except ProcessLookupError:
                pass


@pytest.mark.parametrize("failure", [False, True])
def test_native_sigterm_handler_is_scoped_and_stops_workers(monkeypatch, failure):
    import signal
    original = signal.getsignal(signal.SIGTERM)
    calls = []

    class Process:
        pid = 765432

        def wait(self, timeout=None):
            if timeout is None and failure:
                signal.getsignal(signal.SIGTERM)(signal.SIGTERM, None)
            return 0

    monkeypatch.setattr(earlgrey.subprocess, "Popen", lambda *args, **kwargs: Process())
    monkeypatch.setattr(earlgrey.os, "killpg", lambda pid, sig: calls.append((pid, sig)))
    if failure:
        with pytest.raises(SystemExit) as error:
            earlgrey.run_native(["earlGrey"])
        assert error.value.code == 128 + signal.SIGTERM
        assert (Process.pid, signal.SIGTERM) in calls
        assert (Process.pid, signal.SIGKILL) in calls
    else:
        assert earlgrey.run_native(["earlGrey"]) == 0
        assert not calls
    assert signal.getsignal(signal.SIGTERM) is original
