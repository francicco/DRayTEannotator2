from __future__ import annotations

import subprocess
from dataclasses import dataclass

from .ids import clean_family_id
from pathlib import Path
from typing import Dict, List


@dataclass
class DomainHit:
    family_id: str
    orf_id: str
    domain: str
    accession: str
    evalue: float
    score: float
    hmm_start: int
    hmm_end: int
    ali_start: int
    ali_end: int


def run_hmmscan(
    hmm_db: str | Path,
    proteins_fasta: str | Path,
    domtblout: str | Path,
    hmmscan_bin: str = "hmmscan",
    cpu: int = 1,
) -> Path:
    domtblout = Path(domtblout)
    domtblout.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        hmmscan_bin,
        "--cpu", str(cpu),
        "--domtblout", str(domtblout),
        str(hmm_db),
        str(proteins_fasta),
    ]

    with open(domtblout.with_suffix(".log"), "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)

    return domtblout


def family_from_orf_id(orf_id: str) -> str:
    if "_orf" in orf_id:
        return clean_family_id(orf_id.rsplit("_orf", 1)[0])
    return clean_family_id(orf_id)


def parse_domtblout(
    domtblout: str | Path,
    max_evalue: float = 1e-5,
    min_score: float = 20.0,
) -> List[DomainHit]:
    hits: List[DomainHit] = []

    with open(domtblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split()
            if len(parts) < 23:
                continue

            domain = parts[0]
            accession = parts[1]
            orf_id = parts[3]

            evalue = float(parts[12])
            score = float(parts[13])

            if evalue > max_evalue:
                continue
            if score < min_score:
                continue

            hits.append(
                DomainHit(
                    family_id=family_from_orf_id(orf_id),
                    orf_id=orf_id,
                    domain=domain,
                    accession=accession,
                    evalue=evalue,
                    score=score,
                    hmm_start=int(parts[15]),
                    hmm_end=int(parts[16]),
                    ali_start=int(parts[17]),
                    ali_end=int(parts[18]),
                )
            )

    return hits


def summarize_domains_by_family(hits: List[DomainHit]) -> Dict[str, dict]:
    summary: Dict[str, dict] = {}

    for hit in hits:
        fam = summary.setdefault(
            hit.family_id,
            {
                "domains": set(),
                "best_domain_score": 0.0,
                "domain_hits": 0,
            },
        )

        fam["domains"].add(hit.domain)
        fam["best_domain_score"] = max(fam["best_domain_score"], hit.score)
        fam["domain_hits"] += 1

    for fam in summary.values():
        fam["domains"] = sorted(fam["domains"])

    return summary
