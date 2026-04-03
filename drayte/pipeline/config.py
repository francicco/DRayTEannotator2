from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import yaml


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
    with open(config_file, "r") as handle:
        raw = yaml.safe_load(handle)

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
