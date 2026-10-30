import json
import logging
import signal
import subprocess

import pytest

from drayte.pipeline import reclassify


class InterruptedOutput:
    def __init__(self, error):
        self.error = error
        self.closed = False

    def __iter__(self):
        raise self.error

    def close(self):
        self.closed = True


def run_with_process(tmp_path, monkeypatch, process, killpg):
    merged = tmp_path / "merged.fa"
    merged.write_text(">family#Unknown\nACGT\n")
    calls = []

    def popen(*args, **kwargs):
        calls.append(kwargs)
        return process

    monkeypatch.setattr(reclassify.subprocess, "Popen", popen)
    monkeypatch.setattr(reclassify.os, "killpg", killpg)
    return lambda: reclassify.run_repeatclassifier(
        merged, tmp_path, "RepeatClassifier", logging.getLogger("test.interrupt")
    ), calls


def assert_incomplete(tmp_path):
    record = json.loads((tmp_path / "repeatclassifier.checkpoint.json").read_text())
    assert record["status"] == "running"


def test_ctrl_c_terminates_own_process_group_before_reraising(tmp_path, monkeypatch):
    original = KeyboardInterrupt("stop")
    signals = []

    class Process:
        pid = 424242
        stdout = InterruptedOutput(original)
        waits = []

        def wait(self, timeout=None):
            self.waits.append(timeout)
            return -signal.SIGTERM

    def killpg(pid, sig):
        signals.append((pid, sig))
        if sig == 0:
            raise ProcessLookupError

    process = Process()
    run, calls = run_with_process(tmp_path, monkeypatch, process, killpg)
    with pytest.raises(KeyboardInterrupt) as raised:
        run()
    assert raised.value is original
    assert calls[0]["start_new_session"] is True
    assert signals == [(process.pid, signal.SIGTERM), (process.pid, 0)]
    assert process.waits == [10]
    assert process.stdout.closed
    assert_incomplete(tmp_path)


def test_termination_timeout_kills_group_and_preserves_original_interrupt(tmp_path, monkeypatch):
    original = KeyboardInterrupt("stop")
    signals = []

    class Process:
        pid = 424242
        stdout = InterruptedOutput(original)
        waits = []

        def wait(self, timeout=None):
            self.waits.append(timeout)
            if len(self.waits) == 1:
                raise subprocess.TimeoutExpired("RepeatClassifier", timeout)
            return -signal.SIGKILL

    process = Process()
    run, _ = run_with_process(tmp_path, monkeypatch, process, lambda pid, sig: signals.append((pid, sig)))
    with pytest.raises(KeyboardInterrupt) as raised:
        run()
    assert raised.value is original
    assert signals == [(process.pid, signal.SIGTERM), (process.pid, signal.SIGKILL)]
    assert process.waits == [10, 10]
    assert_incomplete(tmp_path)


def test_nonzero_exit_stops_surviving_descendants(tmp_path, monkeypatch):
    signals = []

    class Process:
        pid = 424242
        stdout = []
        waits = []

        def wait(self, timeout=None):
            self.waits.append(timeout)
            return 7

    process = Process()
    run, _ = run_with_process(tmp_path, monkeypatch, process, lambda pid, sig: signals.append((pid, sig)))
    with pytest.raises(RuntimeError, match="exit code 7"):
        run()
    assert signals == [(process.pid, signal.SIGTERM), (process.pid, 0), (process.pid, signal.SIGKILL)]
    assert process.waits == [None, 10]
    assert_incomplete(tmp_path)


def test_cleanup_failure_does_not_replace_original_exception(tmp_path, monkeypatch):
    original = KeyboardInterrupt("original")

    class Process:
        pid = 424242
        stdout = InterruptedOutput(original)

        def wait(self, timeout=None):
            return -signal.SIGKILL

    def denied(pid, sig):
        raise PermissionError("simulated cleanup failure")

    run, _ = run_with_process(tmp_path, monkeypatch, Process(), denied)
    with pytest.raises(KeyboardInterrupt) as raised:
        run()
    assert raised.value is original
    assert_incomplete(tmp_path)


def test_missing_group_is_harmless_and_wait_reaps_child(tmp_path, monkeypatch):
    original = KeyboardInterrupt("stop")

    class Process:
        pid = 424242
        stdout = InterruptedOutput(original)
        waits = []

        def wait(self, timeout=None):
            self.waits.append(timeout)
            return 1

    def missing(pid, sig):
        raise ProcessLookupError

    process = Process()
    run, _ = run_with_process(tmp_path, monkeypatch, process, missing)
    with pytest.raises(KeyboardInterrupt):
        run()
    assert process.waits == [10]
    assert_incomplete(tmp_path)


def test_output_close_failure_does_not_replace_original_interrupt(tmp_path, monkeypatch):
    original = KeyboardInterrupt("original")

    class Output(InterruptedOutput):
        def close(self):
            raise OSError("simulated close failure")

    class Process:
        pid = 424242
        stdout = Output(original)

        def wait(self, timeout=None):
            return 1

    def missing(pid, sig):
        raise ProcessLookupError

    run, _ = run_with_process(tmp_path, monkeypatch, Process(), missing)
    with pytest.raises(KeyboardInterrupt) as raised:
        run()
    assert raised.value is original
    assert_incomplete(tmp_path)
