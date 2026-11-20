import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from drayte.pipeline import runner
from drayte.pipeline.config import PipelineConfig
from drayte.pipeline.preflight import PreflightError, check_preflight


def write_file(path, text="present\n", *, executable=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    if executable:
        path.chmod(0o755)
    return path


def fake_tool(path):
    # These tools must only be located by preflight, never executed.
    return write_file(path, "#!/bin/sh\nexit 97\n", executable=True)


@pytest.fixture
def ready_config(tmp_path, monkeypatch):
    tools = tmp_path / "bin"
    for name in ("perl", "makeblastdb", "blastn", "faToTwoBit", "RepeatMasker", "RepeatClassifier", "diamond"):
        fake_tool(tools / name)
    monkeypatch.setenv("PATH", str(tools))
    modeler = tmp_path / "RepeatModeler"
    fake_tool(modeler / "util/Linup")
    for name in ("util/alignAndCallConsensus.pl", "util/visualizeAlignPNG.pl", "Matrices/crossmatch/18p41g.matrix", "Matrices/ncbi/nt/25p41g-Hpad.matrix"):
        write_file(modeler / name)
    ucsc = tmp_path / "ucsc"
    for name in ("faToTwoBit", "twoBitInfo", "twoBitToFa"):
        fake_tool(ucsc / name)
    fake_tool(tmp_path / "RepeatAfterMe/RAMExtend")
    fake_tool(tmp_path / "phrap/cross_match")
    write_file(tmp_path / "repeatpeps/RepeatPeps.lib", ">protein\nMAAA\n")
    genome = write_file(tmp_path / "genome.fa", ">chr1\nACGT\n")
    library = write_file(tmp_path / "library.fa", ">fam1#Unknown\nACGT\n")
    script = write_file(tmp_path / "extend.pl", "# invoked with perl; executable permission is unnecessary\n")
    return PipelineConfig(
        genome=str(genome), outdir=str(tmp_path / "results"), species="Mmex",
        extra={
            "repeatmodeler_dir": str(modeler),
            "repeatscout_dir": str(tmp_path / "unused-RepeatScout"),
            "repeatmodeler_library": str(library),
            "repeatmodeler_extend_script": str(script),
            "ucsc_tools_dir": str(ucsc),
            "repeat_after_me_dir": str(tmp_path / "RepeatAfterMe"),
            "phrap_dir": str(tmp_path / "phrap"),
            "repeatpeps_db_dir": str(tmp_path / "repeatpeps"),
        },
    )


def config_yaml(config, tmp_path):
    path = tmp_path / "Mmex.yaml"
    path.write_text(yaml.safe_dump({"genome": config.genome, "outdir": config.outdir,
                                    "species": config.species, **config.extra}))
    return path


def enable_heliano(config, tmp_path):
    for name in ("heliano", "mmseqs", "bedtools", "rnabob", "cd-hit-est", "mafft", "hmmsearch", "getorf", "gt", "dialign2-2", "Rscript", "bash"):
        fake_tool(tmp_path / "bin" / name)
    config.extra["structure"] = {"heliano": {"enabled": True}}


def test_available_workflow_passes_without_running_tools_or_writing_outputs(ready_config):
    report = check_preflight(ready_config)
    assert report.tools and report.inputs
    assert not report.warnings
    assert not ready_config.outdir_path.exists()
    assert not Path(ready_config.genome + ".fai").exists()


def test_reports_all_missing_tools_and_inputs_together(ready_config, tmp_path):
    enable_heliano(ready_config, tmp_path)
    (tmp_path / "bin/heliano").unlink()
    (tmp_path / "bin/Rscript").unlink()
    (tmp_path / "bin/diamond").unlink()
    Path(ready_config.genome).unlink()
    with pytest.raises(PreflightError) as caught:
        check_preflight(ready_config)
    message = str(caught.value)
    for text in ("Genome FASTA", "DIAMOND", "HELIANO", "Rscript", "structure.heliano.bin"):
        assert text in message
    assert not ready_config.outdir_path.exists()


def test_disabled_optional_tools_are_not_required(ready_config):
    ready_config.extra.update({"getorf_bin": "/missing/getorf", "cdhit_est_bin": "/missing/cd-hit-est"})
    ready_config.extra["structure"] = {"heliano": {"enabled": False, "bin": "/missing/heliano"}}
    ready_config.extra["classification"] = {
        "dfam": {"enabled": False, "database": "/missing/dfam"},
        "dfam_db": "/missing/legacy.hmm", "nhmmer_bin": "/missing/nhmmer",
    }
    report = check_preflight(ready_config)
    assert not any("HELIANO" in label or "NHMMER" in label for label, _, _ in report.tools)


def test_enabled_heliano_checks_its_dependencies(ready_config, tmp_path):
    enable_heliano(ready_config, tmp_path)
    (tmp_path / "bin/rnabob").unlink()
    (tmp_path / "bin/getorf").unlink()
    with pytest.raises(PreflightError) as caught:
        check_preflight(ready_config)
    assert "rnabob" in str(caught.value)
    assert "getorf" in str(caught.value)


def test_configured_executable_does_not_fall_back_to_path(ready_config, tmp_path):
    enable_heliano(ready_config, tmp_path)
    missing = str(tmp_path / "wrong-environment/heliano")
    ready_config.extra["structure"]["heliano"]["bin"] = missing
    with pytest.raises(PreflightError) as caught:
        check_preflight(ready_config)
    assert missing in str(caught.value)
    assert "structure.heliano.bin" in str(caught.value)


def test_existing_non_executable_tool_fails(ready_config, tmp_path):
    (tmp_path / "phrap/cross_match").chmod(0o644)
    with pytest.raises(PreflightError, match="cross_match"):
        check_preflight(ready_config)


def test_supplied_library_skips_discovery_tools_but_keeps_extension_helpers(ready_config, tmp_path):
    check_preflight(ready_config)
    (tmp_path / "RepeatModeler/util/Linup").unlink()
    with pytest.raises(PreflightError, match="Linup"):
        check_preflight(ready_config)


def test_discovery_without_supplied_library_requires_modeler_and_scout(ready_config):
    del ready_config.extra["repeatmodeler_library"]
    with pytest.raises(PreflightError) as caught:
        check_preflight(ready_config)
    for name in ("BuildDatabase", "RepeatModeler", "RepeatScout", "build_lmer_table"):
        assert name in str(caught.value)


def test_discovery_tools_can_be_checked_without_execution(ready_config, tmp_path):
    del ready_config.extra["repeatmodeler_library"]
    for name in ("BuildDatabase", "RepeatModeler"):
        fake_tool(tmp_path / "RepeatModeler" / name)
    for name in ("RepeatScout", "build_lmer_table"):
        fake_tool(Path(ready_config.extra["repeatscout_dir"]) / name)
    check_preflight(ready_config)


def test_pfam_requires_pressed_indexes_and_hmmscan(ready_config, tmp_path):
    hmm = write_file(tmp_path / "Pfam-A.hmm")
    ready_config.extra["classification"] = {"pfam_db": str(hmm)}
    with pytest.raises(PreflightError) as caught:
        check_preflight(ready_config)
    for text in ("HMMSCAN", ".h3f", ".h3i", ".h3m", ".h3p", "hmmpress"):
        assert text in str(caught.value)
    fake_tool(tmp_path / "bin/hmmscan")
    for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
        write_file(Path(str(hmm) + suffix))
    check_preflight(ready_config)


def test_famdb_requires_partitions_nhmmer_and_famdb(ready_config, tmp_path):
    data = tmp_path / "FamDB"
    data.mkdir()
    ready_config.extra["classification"] = {"dfam": {"enabled": True, "database": str(data), "taxon": "Formicidae"}}
    with pytest.raises(PreflightError) as caught:
        check_preflight(ready_config)
    for text in ("NHMMER", "FamDB", "*.h5"):
        assert text in str(caught.value)
    for name in ("nhmmer", "famdb.py"):
        fake_tool(tmp_path / "bin" / name)
    write_file(data / "dfam40.0.h5")
    check_preflight(ready_config)


def test_exported_dfam_hmm_does_not_require_famdb_or_pressed_indexes(ready_config, tmp_path):
    hmm = write_file(tmp_path / "selected.hmm")
    fake_tool(tmp_path / "bin/nhmmer")
    ready_config.extra["classification"] = {"dfam": {"enabled": True, "database": str(hmm), "taxon": None}}
    report = check_preflight(ready_config)
    assert not any(label == "FamDB" for label, _, _ in report.tools)


def test_dfam_warn_policy_stays_a_warning(ready_config, tmp_path):
    ready_config.extra["classification"] = {"dfam": {"enabled": True, "database": str(tmp_path / "missing.hmm"), "on_missing": "warn"}}
    report = check_preflight(ready_config)
    assert any("NHMMER" in warning for warning in report.warnings)
    assert any("Dfam HMM" in warning for warning in report.warnings)


@pytest.mark.parametrize("dfam", [{"enabled": True, "engine": "hmmscan"}, {"enabled": True, "threads": 0}])
def test_invalid_dfam_settings_remain_fatal_under_warn(ready_config, dfam):
    ready_config.extra["classification"] = {"dfam": {"on_missing": "warn", **dfam}}
    with pytest.raises(PreflightError, match="classification.dfam"):
        check_preflight(ready_config)


def test_dfam_directory_needs_taxon_and_exported_hmm_disallows_taxon(ready_config, tmp_path):
    data = tmp_path / "FamDB"
    write_file(data / "dfam40.0.h5")
    for name in ("nhmmer", "famdb.py"):
        fake_tool(tmp_path / "bin" / name)
    cfg = {"enabled": True, "database": str(data)}
    ready_config.extra["classification"] = {"dfam": cfg}
    with pytest.raises(PreflightError, match="set the lineage"):
        check_preflight(ready_config)
    cfg.update(database=str(write_file(tmp_path / "selected.hmm")), taxon="Formicidae")
    with pytest.raises(PreflightError, match="use null"):
        check_preflight(ready_config)


def test_dfam_symlink_is_classified_by_resolved_reference(ready_config, tmp_path):
    hmm = write_file(tmp_path / "exported.hmm")
    link = tmp_path / "selected.h5"
    link.symlink_to(hmm)
    fake_tool(tmp_path / "bin/nhmmer")
    ready_config.extra["classification"] = {"dfam": {"enabled": True, "database": str(link)}}
    check_preflight(ready_config)
    partition = write_file(tmp_path / "partition.h5")
    link.unlink()
    link.symlink_to(partition)
    with pytest.raises(PreflightError, match="not an individual .h5"):
        check_preflight(ready_config)


def test_legacy_dfam_and_mmseqs_reference_are_checked(ready_config, tmp_path):
    ready_config.extra["classification"] = {"dfam_db": str(tmp_path / "legacy.hmm"), "mmseqs_rescue_db": str(tmp_path / "rescue.fa")}
    with pytest.raises(PreflightError) as caught:
        check_preflight(ready_config)
    for text in ("NHMMER", "Classification MMseqs", "legacy.hmm", "rescue.fa"):
        assert text in str(caught.value)


def test_missing_repeatpeps_fails_even_with_diamond_database_and_wget(ready_config, tmp_path):
    (tmp_path / "repeatpeps/RepeatPeps.lib").unlink()
    write_file(tmp_path / "repeatpeps/RepeatPeps.lib.dmnd")
    fake_tool(tmp_path / "bin/wget")
    with pytest.raises(PreflightError, match="RepeatPeps.lib"):
        check_preflight(ready_config)
    assert not (tmp_path / "repeatpeps/RepeatPeps.lib").exists()


def test_home_relative_rescue_reference_is_supported(ready_config, tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    write_file(tmp_path / "references/curated.fa", ">ref#DNA\nACGT\n")
    fake_tool(tmp_path / "bin/mmseqs")
    ready_config.extra["classification"] = {"mmseqs_rescue_db": "~/references/curated.fa"}
    report = check_preflight(ready_config)
    assert any(label == "MMseqs rescue FASTA" for label, _, _ in report.inputs)


@pytest.mark.parametrize("database", [42, "empty.hmm"])
def test_invalid_or_empty_dfam_reference_is_not_hidden_by_warn_policy(ready_config, tmp_path, database):
    if database == "empty.hmm":
        database = str(write_file(tmp_path / database, ""))
    fake_tool(tmp_path / "bin/nhmmer")
    ready_config.extra["classification"] = {"dfam": {"enabled": True, "database": database, "on_missing": "warn"}}
    with pytest.raises(PreflightError):
        check_preflight(ready_config)


def test_optional_renderer_and_ucsc_hints_are_advisory(ready_config, tmp_path):
    (tmp_path / "RepeatModeler/util/visualizeAlignPNG.pl").unlink()
    (tmp_path / "ucsc/twoBitInfo").unlink()
    report = check_preflight(ready_config)
    assert any("renderer" in warning for warning in report.warnings)
    assert any("RepModelConfig" in warning for warning in report.warnings)


@pytest.mark.parametrize("key", ["repeatmodeler_dir", "phrap_dir", "repeat_after_me_dir", "repeatmodeler_extend_script", "repeatclassifier_bin"])
def test_paths_used_after_working_directory_changes_must_be_absolute(ready_config, key):
    ready_config.extra[key] = "relative/tool"
    with pytest.raises(PreflightError, match="absolute"):
        check_preflight(ready_config)


def test_relative_output_directory_is_reported_before_execution(ready_config):
    ready_config.outdir = "relative-results"
    with pytest.raises(PreflightError, match="outdir: use an absolute"):
        check_preflight(ready_config)


def test_relative_path_entry_is_not_accepted_for_extension_perl(ready_config, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PATH", "bin")
    with pytest.raises(PreflightError, match="perl: resolves through a relative PATH entry"):
        check_preflight(ready_config)


def test_check_only_returns_before_any_output_mutation(ready_config, tmp_path, monkeypatch, capsys):
    path = config_yaml(ready_config, tmp_path)
    monkeypatch.setattr(sys, "argv", ["drayte", "--config", str(path), "--check-only"])
    for hook in ("output_lock", "ArchiveStore", "run_pipeline", "setup_logger"):
        monkeypatch.setattr(runner, hook, lambda *a, **k: pytest.fail("check-only attempted a mutation"))
    runner.main()
    assert "Pre-flight check passed" in capsys.readouterr().err
    assert not ready_config.outdir_path.exists()


def test_successful_automatic_preflight_continues_to_pipeline(ready_config, tmp_path, monkeypatch, capsys):
    path = config_yaml(ready_config, tmp_path)
    monkeypatch.setattr(sys, "argv", ["drayte", "--config", str(path)])
    seen = []
    monkeypatch.setattr(runner, "run_pipeline", lambda config, filename: seen.append(config.species))
    runner.main()
    assert seen == ["Mmex"]
    assert "Pre-flight check passed" in capsys.readouterr().err
    assert (ready_config.outdir_path / ".storage/run.lock").is_file()


def test_cli_library_override_is_applied_before_preflight(ready_config, tmp_path, monkeypatch):
    supplied = ready_config.extra.pop("repeatmodeler_library")
    path = config_yaml(ready_config, tmp_path)
    monkeypatch.setattr(sys, "argv", ["drayte", "--config", str(path), "--repeatmodeler-library", supplied, "--check-only"])
    runner.main()  # No BuildDatabase/RepeatModeler/RepeatScout exists in this fixture.
    assert not ready_config.outdir_path.exists()


def test_cli_rescue_override_is_applied_before_preflight(ready_config, tmp_path, monkeypatch):
    path = config_yaml(ready_config, tmp_path)
    rescue = write_file(tmp_path / "reference.fa", ">ref#DNA\nACGT\n")
    monkeypatch.setattr(sys, "argv", ["drayte", "--config", str(path), "--mmseqs-rescue-db", str(rescue), "--check-only"])
    with pytest.raises(SystemExit) as caught:
        runner.main()
    assert "classification.mmseqs_bin" in str(caught.value)


def test_failure_leaves_existing_archives_unchanged(ready_config, tmp_path, monkeypatch):
    enable_heliano(ready_config, tmp_path)
    (tmp_path / "bin/heliano").unlink()
    write_file(ready_config.outdir_path / ".storage/archives.json", '{"must": "stay untouched"}\n')
    write_file(ready_config.outdir_path / ".storage/evidence.tar.gz", "archive payload")
    before = {p.relative_to(ready_config.outdir_path): p.read_bytes() for p in ready_config.outdir_path.rglob("*") if p.is_file()}
    path = config_yaml(ready_config, tmp_path)
    monkeypatch.setattr(sys, "argv", ["drayte", "--config", str(path)])
    with pytest.raises(SystemExit, match="Pre-flight check failed"):
        runner.main()
    after = {p.relative_to(ready_config.outdir_path): p.read_bytes() for p in ready_config.outdir_path.rglob("*") if p.is_file()}
    assert after == before


def test_missing_heliano_cli_reports_setting_without_traceback_or_new_output(ready_config, tmp_path):
    enable_heliano(ready_config, tmp_path)
    (tmp_path / "bin/heliano").unlink()
    path = config_yaml(ready_config, tmp_path)
    result = subprocess.run(
        [sys.executable, "-m", "drayte.pipeline.runner", "--config", str(path)],
        cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "Pre-flight check failed" in result.stderr
    assert "structure.heliano.bin" in result.stderr
    assert "Traceback" not in result.stderr
    assert "STAGE:" not in result.stdout + result.stderr
    assert not ready_config.outdir_path.exists()
