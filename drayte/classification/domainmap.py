from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Set

import yaml


DEFAULT_ONTOLOGY = (
    Path(__file__).resolve().parents[2]
    / "config"
    / "domain_ontology.yaml"
)


@lru_cache(maxsize=1)
def load_domain_ontology(path: str | Path | None = None):
    ontology_path = Path(path) if path else DEFAULT_ONTOLOGY

    with open(ontology_path) as fh:
        return yaml.safe_load(fh)


def normalize_domain_name(
    domain: str,
    ontology_path: str | Path | None = None,
) -> str | None:

    ontology = load_domain_ontology(ontology_path)

    for normalized, config in ontology.items():

        for alias in config.get("aliases", []):

            pattern = re.escape(alias)

            if re.search(pattern, domain, flags=re.IGNORECASE):
                return normalized

    return None


def normalize_domains(
    domains: Iterable[str],
    ontology_path: str | Path | None = None,
) -> Set[str]:

    normalized = set()

    for domain in domains:

        value = normalize_domain_name(
            domain,
            ontology_path=ontology_path,
        )

        if value:
            normalized.add(value)

    return normalized
