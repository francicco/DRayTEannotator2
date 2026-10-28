"""Exercise family scheduling with real output collection and no external tools."""

import csv
import json
import logging
import threading
from pathlib import Path
from types import SimpleNamespace

import pytest

from drayte.pipeline import extension as extension_stage
from drayte.utils.names import safe_filename


@pytest.fixture
def extension_run(tmp_path, monkeypatch):
    def prepare(name="run", family_ids=None, threads=2, **extra):
        family_ids = family_ids or ["T-fam-1#Unknown", "T-fam-2#LTR_Gypsy", "T-fam-3#LINE"]
        root = tmp_path / name
        assemblies = root / "discovery" / "assemblies_dir"
        rmodeler = root / "discovery" / "rmodeler_dir"
        assemblies.mkdir(parents=True)
        rmodeler.mkdir(parents=True)
        genome = assemblies / "T.fa"
        genome.write_text(">chr1\n" + "ACGT" * 50 + "\n")
        (assemblies / "T.2bit").write_text("existing two-bit genome\n")
        candidate = rmodeler / "consensi.fa.classified"
        fasta = "".join(f">{te_id}\n{'ACGT' * 30}\n" for te_id in family_ids)
        candidate.write_text(fasta)
        (rmodeler / "T-families.mod.fa").write_text(fasta)
        script = tmp_path / "davidExtendConsRAM.pl"
        script.write_text("#!/usr/bin/env perl\n")
        config = SimpleNamespace(
            outdir_path=root,
            genome_path=genome,
            species="T",
            threads=threads,
            extra={
                "repeatmodeler_extend_script": str(script),
                "ucsc_tools_dir": str(tmp_path / "kent"),
                "repeatmodeler_dir": str(tmp_path / "RepeatModeler"),
                "repeat_after_me_dir": str(tmp_path / "RepeatAfterMe"),
                "phrap_dir": str(tmp_path / "phrap"),
                **extra,
            },
        )
        extracted = {}
        for index, te_id in enumerate(family_ids):
            path = root / "family_inputs" / f"family-{index}.fa"
            path.parent.mkdir(exist_ok=True)
            path.write_text(f">{te_id}\n{'ACGT' * 30}\n")
            extracted[te_id] = path

        calls = []
        failures = set()
        lock = threading.Lock()

        def write_outputs(**kwargs):
            workdir = Path(kwargs["workdir"])
            workdir.mkdir(parents=True, exist_ok=True)
            te_id = Path(kwargs["family_fasta"]).read_text().splitlines()[0][1:]
            with lock:
                calls.append(te_id)
            if te_id in failures:
                raise RuntimeError(f"simulated failure: {te_id}")
            (workdir / "rep").write_text(">repam-newrep\n" + "ACGT" * 30 + "\n")
            (workdir / "repseq.unextended").write_text(
                ">CONSENSUS-original\nACGT\n"
                + "".join(f">hit-{index}\nACGT\n" for index in range(12))
            )
            (workdir / "MSA-extended_with_rmod_cons.fa").write_text(
                ">CORECONS\nACGT\n>repam-newrep\nACGT\n"
            )
            (workdir / "img.png").write_bytes(b"family image")
            Path(kwargs["logfile"]).write_text(f"extended {te_id}\n")

        monkeypatch.setattr(extension_stage, "run_blast", lambda **kwargs: root / "blast.out")
        monkeypatch.setattr(extension_stage, "run_extract_align", lambda **kwargs: extracted)
        monkeypatch.setattr(extension_stage, "run_extend_consensus", write_outputs)

        def run():
            return extension_stage.run(
                config,
                {"raw_library": str(candidate)},
                logging.getLogger("test.extension.parallel"),
            )

        return SimpleNamespace(
            root=root,
            config=config,
            extracted=extracted,
            calls=calls,
            failures=failures,
            write_outputs=write_outputs,
            run=run,
        )

    return prepare


def read_rows(result):
    with open(result["summary_tsv"]) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.parametrize("setting", ["omitted", None, 2])
def test_families_overlap_and_worker_limit_is_respected(extension_run, monkeypatch, setting):
    options = {} if setting == "omitted" else {"extension_workers": setting}
    families = [f"T-fam-{index}#Unknown" for index in range(6)]
    run = extension_run(family_ids=families, threads=2, **options)
    barrier = threading.Barrier(2)
    lock = threading.Lock()
    active = 0
    peak = 0

    def extend(**kwargs):
        nonlocal active, peak
        with lock:
            active += 1
            peak = max(peak, active)
        try:
            # A serial implementation cannot pass this barrier; six families
            # exercise multiple waves rather than merely filling two slots.
            barrier.wait(timeout=5)
            run.write_outputs(**kwargs)
        finally:
            with lock:
                active -= 1

    monkeypatch.setattr(extension_stage, "run_extend_consensus", extend)
    result = run.run()

    assert result["n_ok"] == 6
    assert result["n_failed"] == 0
    assert peak == 2
    assert [row["te_id"] for row in read_rows(result)] == families
    for te_id in families:
        safe_id = safe_filename(te_id)
        assert (run.root / "extension" / "extensionwork" / safe_id / "rep").is_file()
        assert (run.root / "extension" / "extendlogs" / f"{safe_id}.extend.log").read_text() == f"extended {te_id}\n"


def test_parallel_results_match_serial_even_when_completion_order_differs(extension_run, monkeypatch):
    families = ["T-first#Unknown", "T-second#Unknown", "T-third#Unknown"]
    serial = extension_run(name="serial", family_ids=families, extension_workers=1)
    serial_result = serial.run()
    parallel = extension_run(name="parallel", family_ids=families, extension_workers=2)
    second_completed = threading.Event()
    completed = []
    lock = threading.Lock()

    def extend(**kwargs):
        te_id = Path(kwargs["family_fasta"]).read_text().splitlines()[0][1:]
        if te_id == families[0]:
            assert second_completed.wait(timeout=5)
        parallel.write_outputs(**kwargs)
        with lock:
            completed.append(te_id)
        if te_id == families[1]:
            second_completed.set()

    monkeypatch.setattr(extension_stage, "run_extend_consensus", extend)
    summary_threads = []
    original_write_summary = extension_stage.write_summary

    def write_summary(*args, **kwargs):
        summary_threads.append(threading.current_thread())
        return original_write_summary(*args, **kwargs)

    monkeypatch.setattr(extension_stage, "write_summary", write_summary)
    parallel_result = parallel.run()

    assert completed.index(families[1]) < completed.index(families[0])
    assert summary_threads == [threading.current_thread()]
    serial_text = Path(serial_result["summary_tsv"]).read_text().replace(str(serial.root), "ROOT")
    parallel_text = Path(parallel_result["summary_tsv"]).read_text().replace(str(parallel.root), "ROOT")
    assert parallel_text == serial_text
    for path in (serial.root / "extension" / "final_consensuses").iterdir():
        assert path.read_bytes() == (parallel.root / "extension" / "final_consensuses" / path.name).read_bytes()


def test_one_worker_executes_serially(extension_run, monkeypatch):
    run = extension_run(threads=4, extension_workers=1)
    active = False

    def extend(**kwargs):
        nonlocal active
        assert not active
        active = True
        try:
            run.write_outputs(**kwargs)
        finally:
            active = False

    monkeypatch.setattr(extension_stage, "run_extend_consensus", extend)
    result = run.run()
    assert result["n_ok"] == len(run.extracted)
    assert run.calls == list(run.extracted)


def test_failed_family_does_not_stop_others_and_resume_retries_only_failure(extension_run):
    run = extension_run(extension_workers=1)
    families = list(run.extracted)
    run.failures.add(families[1])
    result = run.run()
    rows = read_rows(result)
    assert result["n_ok"] == 2
    assert result["n_failed"] == 1
    assert [row["status"] for row in rows] == ["ok", "failed_extension", "ok"]
    assert rows[1]["error_message"] == f"simulated failure: {families[1]}"

    # The runner normally persists this manifest after the stage returns.
    manifest = run.root / "extension" / "extension.manifest.json"
    manifest.write_text(json.dumps(result))
    run.calls.clear()
    run.failures.clear()
    run.config.extra["extension_workers"] = 3
    resumed = run.run()

    assert run.calls == [families[1]]
    assert resumed["n_ok"] == 3
    assert resumed["n_failed"] == 0
    assert resumed["input_signature"] == result["input_signature"]
    assert [row["te_id"] for row in read_rows(resumed)] == families


def test_completed_checkpoint_is_reused_after_changing_worker_count(extension_run, monkeypatch):
    run = extension_run(extension_workers=2)
    result = run.run()
    manifest = run.root / "extension" / "extension.manifest.json"
    manifest.write_text(json.dumps(result))
    run.calls.clear()
    run.config.extra["extension_workers"] = 1

    def unexpected(**kwargs):
        pytest.fail("completed checkpoint must skip BLAST and family extraction")

    monkeypatch.setattr(extension_stage, "run_blast", unexpected)
    monkeypatch.setattr(extension_stage, "run_extract_align", unexpected)
    resumed = run.run()

    assert run.calls == []
    assert resumed["resume_validation"]["input_signature_status"] == "matching_signature"


def test_all_failures_write_resumable_manifest_and_raise(extension_run):
    run = extension_run(extension_workers=2)
    run.failures.update(run.extracted)
    with pytest.raises(RuntimeError, match="no successful consensuses"):
        run.run()
    manifest = run.root / "extension" / "extension.manifest.json"
    failed = json.loads(manifest.read_text())
    assert failed["n_ok"] == 0
    assert failed["n_failed"] == len(run.extracted)
    assert len(read_rows(failed)) == len(run.extracted)

    run.failures.clear()
    resumed = run.run()
    assert resumed["n_ok"] == len(run.extracted)


@pytest.mark.parametrize("workers", [0, -1, True, False, 1.5, "2", [], {}])
def test_invalid_worker_count_is_rejected(extension_run, workers):
    run = extension_run(extension_workers=workers)
    with pytest.raises(ValueError, match="extension_workers"):
        run.run()
    assert run.calls == []


def test_safe_filename_collision_is_rejected_before_jobs_launch(extension_run, monkeypatch):
    run = extension_run(family_ids=["T-family/a", "T-family_a"], extension_workers=2)

    def unexpected(**kwargs):
        pytest.fail("colliding family IDs must be rejected before extraction")

    monkeypatch.setattr(extension_stage, "run_extract_align", unexpected)
    with pytest.raises(ValueError, match="same extension filename"):
        run.run()
    assert run.calls == []
