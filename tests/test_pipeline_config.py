import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from drayte.pipeline.config import ConfigError, load_config


def valid_config(tmp_path):
    return {
        "genome": str(tmp_path / "genome.fa"),
        "outdir": str(tmp_path / "results"),
        "species": "Mmex",
    }


@pytest.mark.parametrize("key", ["genome", "outdir", "species"])
@pytest.mark.parametrize("value", [None, "", "   ", 42, False, [], {}])
def test_required_fields_report_invalid_values(tmp_path, key, value):
    raw = valid_config(tmp_path)
    raw[key] = value
    path = tmp_path / "Mmex.yaml"
    path.write_text(yaml.safe_dump(raw))

    with pytest.raises(ConfigError) as caught:
        load_config(str(path))

    message = str(caught.value)
    assert str(path) in message
    assert f"'{key}'" in message
    assert "Set" in message
    if key == "genome":
        assert "genome FASTA" in message
        assert "GENOME" in message
    assert not (tmp_path / "results").exists()


@pytest.mark.parametrize("key", ["genome", "outdir", "species"])
def test_missing_required_field_names_key_and_config_file(tmp_path, key):
    raw = valid_config(tmp_path)
    del raw[key]
    path = tmp_path / "Mmex.yaml"
    path.write_text(yaml.safe_dump(raw))

    with pytest.raises(ConfigError) as caught:
        load_config(str(path))

    assert str(path) in str(caught.value)
    assert f"'{key}' is missing or empty" in str(caught.value)


@pytest.mark.parametrize("text", ["", "# empty config\n", "- genome.fa\n", "42\n", "true\n"])
def test_non_mapping_yaml_reports_expected_configuration(tmp_path, text):
    path = tmp_path / "Mmex.yaml"
    path.write_text(text)

    with pytest.raises(ConfigError, match="expected a YAML mapping") as caught:
        load_config(str(path))

    assert str(path) in str(caught.value)


def test_malformed_yaml_reports_file_and_location(tmp_path):
    path = tmp_path / "Mmex.yaml"
    path.write_text("genome: [\noutdir: results\n")

    with pytest.raises(ConfigError) as caught:
        load_config(str(path))

    message = str(caught.value)
    assert str(path) in message
    assert "Invalid YAML" in message
    assert "line " in message and "column " in message
    assert "\n" not in message


@pytest.mark.parametrize("source", ["missing", "directory", "invalid_utf8"])
def test_unreadable_config_reports_filename(tmp_path, source):
    path = tmp_path / "Mmex.yaml"
    if source == "directory":
        path.mkdir()
    elif source == "invalid_utf8":
        path.write_bytes(b"\xff\xfe")

    with pytest.raises(ConfigError, match="Cannot read configuration file") as caught:
        load_config(str(path))

    assert str(path) in str(caught.value)


def test_valid_configuration_preserves_values_and_extra_settings(tmp_path):
    raw = valid_config(tmp_path)
    raw["genome"] = " genome with spaces.fa "
    raw["threads"] = 16
    raw["storage"] = {"mode": "compact"}
    raw["classification"] = {"dfam": {"enabled": True, "taxon": "Formicidae"}}
    path = tmp_path / "Mmex.yaml"
    path.write_text(yaml.safe_dump(raw))

    config = load_config(str(path))

    assert config.genome == raw["genome"]
    assert config.outdir == raw["outdir"]
    assert config.species == "Mmex"
    assert config.threads == 16
    assert config.extra == {"storage": raw["storage"], "classification": raw["classification"]}


def test_empty_genome_cli_exits_without_traceback_or_output_directory(tmp_path):
    raw = valid_config(tmp_path)
    raw["genome"] = None
    path = tmp_path / "Mmex.drayannConfig.Assisted.yaml"
    path.write_text(yaml.safe_dump(raw))

    result = subprocess.run(
        [sys.executable, "-m", "drayte.pipeline.runner", "--config", str(path)],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert result.stderr.startswith("drayte: error: ")
    assert str(path) in result.stderr
    assert "'genome' is missing or empty" in result.stderr
    assert "full path of the genome FASTA file" in result.stderr
    assert "GENOME" in result.stderr
    assert "Traceback" not in result.stderr
    assert "STAGE:" not in result.stdout + result.stderr
    assert not Path(raw["outdir"]).exists()
