from __future__ import annotations

from pathlib import Path
from typing import Iterator

from Bio import SeqIO

from .ids import clean_family_id
from .tsd import (
    TSDHit,
    detect_tsd_for_copy,
    summarize_tsd_hits,
)


def extract_flanks(
    seq: str,
    flank_size: int = 15,
) -> tuple[str, str]:
    seq = seq.upper()

    if len(seq) < flank_size * 2:
        return "", ""

    left = seq[:flank_size]
    right = seq[-flank_size:]

    return left, right


def iter_alignment_sequences(
    fasta: str | Path,
) -> Iterator[tuple[str, str]]:
    for rec in SeqIO.parse(str(fasta), "fasta"):
        yield rec.id, str(rec.seq)


def detect_family_tsd_from_alignment(
    family_id: str,
    alignment_fasta: str | Path,
    flank_size: int = 15,
    min_len: int = 2,
    max_len: int = 15,
    min_identity: float = 1.0,
    min_support_fraction: float = 0.5,
) -> TSDHit:
    hits: list[TSDHit] = []

    for _, seq in iter_alignment_sequences(alignment_fasta):

        left, right = extract_flanks(
            seq,
            flank_size=flank_size,
        )

        if not left or not right:
            continue

        hit = detect_tsd_for_copy(
            family_id=family_id,
            left_flank=left,
            right_flank=right,
            min_len=min_len,
            max_len=max_len,
            min_identity=min_identity,
        )

        hits.append(hit)

    return summarize_tsd_hits(
        family_id=family_id,
        hits=hits,
        min_support_fraction=min_support_fraction,
    )
