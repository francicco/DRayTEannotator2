from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import yaml


class ConfigError(ValueError):
    """A configuration problem that can be explained before starting stages."""


@dataclass
class PipelineConfig:
    genome: str
    outdir: str
    species: str
    threads: int = 8
    preset: str = "generic_insect"
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def genome_path(self) -> Path:
        return Path(self.genome)

    @property
    def outdir_path(self) -> Path:
        return Path(self.outdir)


def load_config(config_file: str) -> PipelineConfig:
    try:
        with open(config_file, "r", encoding="utf-8") as handle:
            raw = yaml.safe_load(handle)
    except (OSError, UnicodeError) as error:
        raise ConfigError(f"Cannot read configuration file '{config_file}': {error}") from error
    except yaml.YAMLError as error:
        mark = getattr(error, "problem_mark", None)
        location = f" at line {mark.line + 1}, column {mark.column + 1}" if mark else ""
        problem = getattr(error, "problem", None) or "invalid YAML syntax"
        raise ConfigError(f"Invalid YAML in '{config_file}'{location}: {problem}") from error

    if not isinstance(raw, dict):
        raise ConfigError(
            f"Invalid configuration in '{config_file}': expected a YAML mapping "
            "containing 'genome', 'outdir', and 'species'."
        )

    guidance = {
        "genome": (
            "Set 'genome' to the full path of the genome FASTA file. "
            "If a shell script generated this YAML, check that GENOME is set."
        ),
        "outdir": "Set 'outdir' to the output directory for this run.",
        "species": "Set 'species' to a species tag, such as Mmex.",
    }
    for key, hint in guidance.items():
        value = raw.get(key)
        if value is None or (isinstance(value, str) and not value.strip()):
            raise ConfigError(
                f"Invalid configuration in '{config_file}': '{key}' is missing or empty. {hint}"
            )
        if not isinstance(value, str):
            raise ConfigError(
                f"Invalid configuration in '{config_file}': '{key}' must be a "
                f"non-empty string (got {type(value).__name__}). {hint}"
            )

    known = {"genome", "outdir", "species", "threads", "preset"}
    extra = {k: v for k, v in raw.items() if k not in known}

    return PipelineConfig(
        genome=raw["genome"],
        outdir=raw["outdir"],
        species=raw["species"],
        threads=raw.get("threads", 8),
        preset=raw.get("preset", "generic_insect"),
        extra=extra,
    )
