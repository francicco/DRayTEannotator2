from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path


@dataclass
class TSDHit:
    family_id: str
    tsd_present: bool
    tsd_seq: str
    tsd_len: int
    identity: float
    n_copies: int = 1


def best_tsd_from_boundaries(
    left_flank: str,
    right_flank: str,
    min_len: int = 2,
    max_len: int = 15,
    min_identity: float = 1.0,
) -> tuple[str, float] | None:
    left_flank = left_flank.upper()
    right_flank = right_flank.upper()

    best = None

    for k in range(max_len, min_len - 1, -1):
        if len(left_flank) < k or len(right_flank) < k:
            continue

        left = left_flank[-k:]
        right = right_flank[:k]

        matches = sum(
            1 for a, b in zip(left, right)
            if a == b
        )

        identity = matches / k

        if identity >= min_identity:
            best = (left, identity)
            break

    return best


def detect_tsd_for_copy(
    family_id: str,
    left_flank: str,
    right_flank: str,
    min_len: int = 2,
    max_len: int = 15,
    min_identity: float = 1.0,
) -> TSDHit:
    best = best_tsd_from_boundaries(
        left_flank=left_flank,
        right_flank=right_flank,
        min_len=min_len,
        max_len=max_len,
        min_identity=min_identity,
    )

    if best is None:
        return TSDHit(
            family_id=family_id,
            tsd_present=False,
            tsd_seq="",
            tsd_len=0,
            identity=0.0,
        )

    seq, identity = best

    return TSDHit(
        family_id=family_id,
        tsd_present=True,
        tsd_seq=seq,
        tsd_len=len(seq),
        identity=identity,
    )


def summarize_tsd_hits(
    family_id: str,
    hits: list[TSDHit],
    min_support_fraction: float = 0.5,
) -> TSDHit:
    positive = [h for h in hits if h.tsd_present]

    if not hits or not positive:
        return TSDHit(
            family_id=family_id,
            tsd_present=False,
            tsd_seq="",
            tsd_len=0,
            identity=0.0,
            n_copies=len(hits),
        )

    seq_counts = Counter(h.tsd_seq for h in positive)
    best_seq, best_n = seq_counts.most_common(1)[0]

    support_fraction = best_n / len(hits)

    if support_fraction < min_support_fraction:
        return TSDHit(
            family_id=family_id,
            tsd_present=False,
            tsd_seq=best_seq,
            tsd_len=len(best_seq),
            identity=support_fraction,
            n_copies=len(hits),
        )

    return TSDHit(
        family_id=family_id,
        tsd_present=True,
        tsd_seq=best_seq,
        tsd_len=len(best_seq),
        identity=support_fraction,
        n_copies=len(hits),
    )
