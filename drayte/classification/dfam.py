from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List

from .ids import clean_family_id


@dataclass
class DfamHit:
    family_id: str
    model_name: str
    accession: str
    evalue: float
    score: float
    ali_start: int
    ali_end: int


def run_nhmmer(
    dfam_db: str | Path,
    consensus_fasta: str | Path,
    tblout: str | Path,
    nhmmer_bin: str = "nhmmer",
    cpu: int = 1,
) -> Path:
    tblout = Path(tblout)
    tblout.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        nhmmer_bin,
        "--cpu", str(cpu),
        "--tblout", str(tblout),
        str(dfam_db),
        str(consensus_fasta),
    ]

    with open(tblout.with_suffix(".log"), "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)

    return tblout


def parse_nhmmer_tblout(
    tblout: str | Path,
    max_evalue: float = 1e-5,
    min_score: float = 20.0,
) -> List[DfamHit]:
    hits: List[DfamHit] = []

    with open(tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split()
            if len(parts) < 14:
                continue

            family_id = clean_family_id(parts[0])
            model_name = parts[2]
            accession = parts[3]
            ali_start = int(parts[6])
            ali_end = int(parts[7])
            evalue = float(parts[12])
            score = float(parts[13])

            if evalue > max_evalue:
                continue

            if score < min_score:
                continue

            hits.append(
                DfamHit(
                    family_id=family_id,
                    model_name=model_name,
                    accession=accession,
                    evalue=evalue,
                    score=score,
                    ali_start=ali_start,
                    ali_end=ali_end,
                )
            )

    return hits
