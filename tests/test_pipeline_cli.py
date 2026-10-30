import logging
import sys
from pathlib import Path

import pytest
import yaml

from drayte.pipeline import runner


class StopBeforeDiscovery(Exception):
    """Prevent the command-line test from launching external tools."""


def capture_pipeline_config(tmp_path, monkeypatch, classification, extra_args):
    config_path = tmp_path / "config.yaml"
    raw = {
        "genome": str(tmp_path / "genome.fa"),
        "outdir": str(tmp_path / "results"),
        "species": "Hmel",
    }
    if classification != "omitted":
        raw["classification"] = classification
    config_path.write_text(yaml.safe_dump(raw))
    monkeypatch.setattr(sys, "argv", ["drayte", "--config", str(config_path), *extra_args])
    monkeypatch.setattr(runner, "setup_logger", lambda **kwargs: logging.getLogger("test.pipeline.cli"))
    captured = {}

    def stop_discovery(config, logger):
        captured["config"] = config
        raise StopBeforeDiscovery

    monkeypatch.setattr(runner.discovery, "run", stop_discovery)
    with pytest.raises(StopBeforeDiscovery):
        runner.main()
    return captured["config"]


def test_rescue_flag_overrides_yaml_and_preserves_other_settings(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config = capture_pipeline_config(
        tmp_path, monkeypatch,
        {"mmseqs_rescue_db": "old.fa", "pfam_db": "Pfam-A.hmm", "cpu": 2},
        ["--mmseqs-rescue-db", "Heliconius.classification.fa"],
    )
    assert config.extra["classification"] == {
        "mmseqs_rescue_db": str(tmp_path / "Heliconius.classification.fa"),
        "pfam_db": "Pfam-A.hmm",
        "cpu": 2,
    }


def test_omitted_rescue_flag_preserves_yaml(tmp_path, monkeypatch):
    settings = {"mmseqs_rescue_db": "configured.fa", "dfam_db": "Dfam.hmm"}
    config = capture_pipeline_config(tmp_path, monkeypatch, settings, [])
    assert config.extra["classification"] == settings


@pytest.mark.parametrize("classification", ["omitted", None])
def test_rescue_flag_initializes_missing_or_null_section(tmp_path, monkeypatch, classification):
    reference = tmp_path / "Heliconius.classification.fa"
    config = capture_pipeline_config(
        tmp_path, monkeypatch, classification,
        ["--mmseqs-rescue-db", str(reference)],
    )
    assert config.extra["classification"] == {"mmseqs_rescue_db": str(reference)}


def test_rescue_flag_expands_home_and_coexists_with_discovery_override(tmp_path, monkeypatch):
    config = capture_pipeline_config(
        tmp_path, monkeypatch, {},
        ["--mmseqs-rescue-db", "~/Heliconius.classification.fa",
         "--repeatmodeler-library", str(tmp_path / "consensi.fa.classified")],
    )
    assert config.extra["classification"]["mmseqs_rescue_db"] == str(
        (Path.home() / "Heliconius.classification.fa").resolve()
    )
    assert config.extra["repeatmodeler_library"] == str(tmp_path / "consensi.fa.classified")


def test_rescue_flag_rejects_non_mapping_section_before_running_tools(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="classification config section must be a mapping"):
        capture_pipeline_config(tmp_path, monkeypatch, "not-a-mapping",
                                ["--mmseqs-rescue-db", "reference.fa"])
    assert not (tmp_path / "results").exists()


def test_storage_flag_is_forwarded_without_changing_classification(tmp_path, monkeypatch):
    settings = {"mmseqs_rescue_db": "configured.fa", "cpu": 2}
    config = capture_pipeline_config(tmp_path, monkeypatch, settings, ["--storage-mode", "compact"])
    assert config.extra["storage"]["mode"] == "compact"
    assert config.extra["classification"] == settings
