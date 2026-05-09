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
    family_id: str = "unknown",
    window: int = 120,
    min_len: int = 10,
    min_identity: float = 0.75,
    logger=None,
) -> TirDetection | None:
    """
    Search for the best terminal inverted repeat using a sliding comparison
    between the 5' terminal window and the reverse complement of the 3'
    terminal window.
    """

    seq = seq.upper()
    n = len(seq)

    if logger:
        logger.info(
            "[TIR] scanning %s len=%d window=%d min_len=%d min_identity=%.2f",
            family_id,
            n,
            window,
            min_len,
            min_identity,
        )

    if n < min_len * 2:
        if logger:
            logger.info("[TIR] no_TIR %s reason=too_short", family_id)
        return None

    terminal_window = min(window, n)

    left = seq[:terminal_window]
    right = str(
        Seq(seq[-terminal_window:]).reverse_complement()
    )

    best = None

    for i in range(0, len(left) - min_len + 1):
        for j in range(0, len(right) - min_len + 1):

            max_k = min(
                len(left) - i,
                len(right) - j,
            )

            for k in range(min_len, max_k + 1):
                lsub = left[i:i + k]
                rsub = right[j:j + k]

                matches = sum(
                    1 for a, b in zip(lsub, rsub)
                    if a == b
                )

                identity = matches / k

                if identity < min_identity:
                    continue

                if best is None:
                    best = (k, identity, i, j)
                    continue

                best_len, best_identity, _, _ = best

                if (
                    k > best_len
                    or (
                        k == best_len
                        and identity > best_identity
                    )
                ):
                    best = (k, identity, i, j)

    if best is None:
        if logger:
            logger.info("[TIR] no_TIR %s", family_id)
        return None

    tir_len, tir_identity, left_offset, right_offset = best

    left_start = left_offset + 1

    right_start = (
        n
        - terminal_window
        + right_offset
        + 1
    )

    if logger:
        logger.info(
            "[TIR] detected %s len=%d identity=%.3f left=%d right=%d",
            family_id,
            tir_len,
            tir_identity,
            left_start,
            right_start,
        )

    return TirDetection(
        family_id=family_id,
        tir_present=True,
        tir_len=tir_len,
        tir_identity=tir_identity,
        left_start=left_start,
        right_start=right_start,
    )


def detect_tirs_from_fasta(
    fasta: str | Path,
    window: int = 120,
    min_len: int = 10,
    min_identity: float = 0.75,
    logger=None,
) -> list[TirDetection]:

    detections: list[TirDetection] = []

    records = list(SeqIO.parse(str(fasta), "fasta"))
    total = len(records)
    
    for i, rec in enumerate(records, start=1):
        family_id = clean_family_id(rec.id)

        det = best_terminal_inverted_repeat(
            str(rec.seq),
            family_id=family_id,
            window=window,
            min_len=min_len,
            min_identity=min_identity,
            logger=None,
        )

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
            detections.append(det)

        if logger and (i == 1 or i % 25 == 0 or i == total):
            logger.info(
                "[TIR] progress %d/%d families scanned; detected=%d",
                i,
                total,
                sum(1 for d in detections if d.tir_present),
            )

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
                "source": (
                    f"heuristic_tir:"
                    f"len={det.tir_len};"
                    f"left={det.left_start};"
                    f"right={det.right_start}"
                ),
            })
