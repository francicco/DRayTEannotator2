from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from .ids import clean_family_id


@dataclass
class TirDetection:
    family_id: str
    tir_present: bool
    tir_len: int
    tir_identity: float
    left_start: int
    right_start: int


def best_terminal_inverted_repeat(
    seq: str,
    window: int = 250,
    min_len: int = 15,
    min_identity: float = 0.80,
) -> TirDetection | None:
    seq = seq.upper()
    n = len(seq)

    if n < min_len * 2:
        return None

    left = seq[: min(window, n)]
    right = str(Seq(seq[-min(window, n):]).reverse_complement())

    best = None

    for k in range(min_len, min(len(left), len(right)) + 1):
        lsub = left[:k]
        rsub = right[:k]

        matches = sum(1 for a, b in zip(lsub, rsub) if a == b)
        identity = matches / k

        if identity >= min_identity:
            best = (k, identity)

    if best is None:
        return None

    tir_len, tir_identity = best

    return TirDetection(
        family_id="",
        tir_present=True,
        tir_len=tir_len,
        tir_identity=tir_identity,
        left_start=1,
        right_start=n - tir_len + 1,
    )


def detect_tirs_from_fasta(
    fasta: str | Path,
    window: int = 250,
    min_len: int = 15,
    min_identity: float = 0.80,
) -> list[TirDetection]:
    detections = []

    for rec in SeqIO.parse(str(fasta), "fasta"):
        det = best_terminal_inverted_repeat(
            str(rec.seq),
            window=window,
            min_len=min_len,
            min_identity=min_identity,
        )

        family_id = clean_family_id(rec.id)

        if det is None:
            detections.append(
                TirDetection(
                    family_id=family_id,
                    tir_present=False,
                    tir_len=0,
                    tir_identity=0.0,
                    left_start=0,
                    right_start=0,
                )
            )
        else:
            det.family_id = family_id
            detections.append(det)

    return detections


def write_tir_structure_tsv(
    detections: list[TirDetection],
    outpath: str | Path,
) -> None:
    fields = [
        "family_id",
        "evidence_type",
        "score",
        "source",
    ]

    with open(outpath, "w") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()

        for det in detections:
            if not det.tir_present:
                continue

            writer.writerow({
                "family_id": det.family_id,
                "evidence_type": "TIR",
                "score": round(det.tir_identity, 3),
                "source": f"heuristic_tir:len={det.tir_len}",
            })
