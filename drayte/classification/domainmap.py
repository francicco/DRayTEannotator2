from __future__ import annotations

import re
from typing import Iterable, Set


DOMAIN_PATTERNS = {
    "RT": [
        r"^RT$",
        r"RT_domain",
        r"Reverse[_ -]?transcriptase",
        r"RVT",
        r"retropepsin",
    ],
    "INTEGRASE": [
        r"^INT$",
        r"INT_domain",
        r"Integrase",
        r"rve",
    ],
    "RNASEH": [
        r"RNase[_ -]?H",
        r"RNaseH",
        r"RH_domain",
    ],
    "GAG": [
        r"^GAG$",
        r"^Gag$",
        r"gag",
    ],
    "ENV": [
        r"^ENV$",
        r"^Env$",
        r"envelope",
    ],
    "TRANSPOSASE": [
        r"Transposase",
        r"DDE",
        r"DDE_Tnp",
        r"hAT",
        r"Tc1",
        r"Mutator",
        r"Mariner",
        r"piggyBac",
        r"PIF",
    ],
    "HELITRON_REP": [
        r"Helitron",
        r"RepHel",
        r"Replicase",
        r"Rolling[_ -]?circle",
    ],
    "HELICASE": [
        r"Helicase",
        r"DEXDc",
        r"HELICc",
    ],
}


def normalize_domain_name(domain: str) -> str | None:
    for normalized, patterns in DOMAIN_PATTERNS.items():
        for pattern in patterns:
            if re.search(pattern, domain, flags=re.IGNORECASE):
                return normalized
    return None


def normalize_domains(domains: Iterable[str]) -> Set[str]:
    normalized = set()

    for domain in domains:
        value = normalize_domain_name(domain)
        if value:
            normalized.add(value)

    return normalized
