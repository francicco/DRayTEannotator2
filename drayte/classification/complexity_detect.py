from __future__ import annotations

from dataclasses import dataclass
from math import log2
from pathlib import Path
from collections import Counter

from Bio import SeqIO

from drayte.classification.ids import clean_family_id


@dataclass
class ComplexityDetection:
    family_id: str
    consensus_len: int = 0

    seq_entropy: float = 0.0
    dominant_1mer_fraction: float = 0.0
    dominant_2mer_fraction: float = 0.0
    dominant_3mer_fraction: float = 0.0
    dominant_4mer_fraction: float = 0.0
    max_homopolymer_run: int = 0

    min_window_entropy: float = 0.0
    max_window_dom1_fraction: float = 0.0
    max_window_homopolymer_run: int = 0

    low_complexity_candidate: bool = False
    simple_repeat_candidate: bool = False
    complexity_class: str = "none"
    complexity_confidence: str = "LOW"


def _clean_seq(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in {"A", "C", "G", "T"})


def shannon_entropy(seq: str) -> float:
    seq = _clean_seq(seq)
    if not seq:
        return 0.0

    n = len(seq)
    counts = Counter(seq)
    entropy = 0.0

    for count in counts.values():
        p = count / n
        entropy -= p * log2(p)

    return round(entropy, 3)


def dominant_1mer_fraction(seq: str) -> float:
    seq = _clean_seq(seq)
    if not seq:
        return 0.0
    return round(max(Counter(seq).values()) / len(seq), 3)


def max_homopolymer_run(seq: str) -> int:
    seq = _clean_seq(seq)
    if not seq:
        return 0

    best = 1
    current = 1

    for a, b in zip(seq, seq[1:]):
        if a == b:
            current += 1
            best = max(best, current)
        else:
            current = 1

    return best


def dominant_periodic_kmer_fraction(seq: str, k: int) -> float:
    seq = _clean_seq(seq)
    if len(seq) < k * 3:
        return 0.0

    best = 0.0

    for offset in range(k):
        kmers = [
            seq[i:i + k]
            for i in range(offset, len(seq) - k + 1, k)
            if len(seq[i:i + k]) == k
        ]

        if not kmers:
            continue

        counts = Counter(kmers)
        dominant = max(counts.values())
        fraction = dominant / len(kmers)
        best = max(best, fraction)

    return round(best, 3)


def detect_complexity_structure(
    family_id: str,
    seq: str,
) -> ComplexityDetection:
    clean = _clean_seq(seq)
    consensus_len = len(clean)

    entropy = shannon_entropy(clean)
    d1 = dominant_1mer_fraction(clean)
    d2 = dominant_periodic_kmer_fraction(clean, 2)
    d3 = dominant_periodic_kmer_fraction(clean, 3)
    d4 = dominant_periodic_kmer_fraction(clean, 4)
    homo = max_homopolymer_run(clean)
    min_win_entropy, max_win_dom1, max_win_homo = local_complexity_metrics(clean)
    max_periodic = max(d2, d3, d4)

    simple_repeat = bool(
        consensus_len >= 30
        and (
            max_periodic >= 0.70
            or d1 >= 0.80
            or homo >= 20
        )
    )

    low_complexity = bool(
        consensus_len >= 30
        and not simple_repeat
        and (
            entropy <= 1.35
            or d1 >= 0.70
            or max_periodic >= 0.55
            or min_win_entropy <= 1.00
            or max_win_dom1 >= 0.80
            or max_win_homo >= 18
        )
    )

    if simple_repeat:
        complexity_class = "simple_repeat"
    elif low_complexity:
        complexity_class = "low_complexity"
    else:
        complexity_class = "none"

    if simple_repeat and (max_periodic >= 0.85 or d1 >= 0.90 or homo >= 30):
        confidence = "HIGH"
    elif simple_repeat or low_complexity:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"

    return ComplexityDetection(
        family_id=clean_family_id(family_id),
        consensus_len=consensus_len,
        seq_entropy=entropy,
        dominant_1mer_fraction=d1,
        dominant_2mer_fraction=d2,
        dominant_3mer_fraction=d3,
        dominant_4mer_fraction=d4,
        max_homopolymer_run=homo,
        low_complexity_candidate=low_complexity,
        min_window_entropy=min_win_entropy,
        max_window_dom1_fraction=max_win_dom1,
        max_window_homopolymer_run=max_win_homo,
        simple_repeat_candidate=simple_repeat,
        complexity_class=complexity_class,
        complexity_confidence=confidence,
    )


def detect_complexity_from_fasta(
    fasta: str | Path,
) -> list[ComplexityDetection]:
    calls: list[ComplexityDetection] = []

    for record in SeqIO.parse(str(fasta), "fasta"):
        calls.append(
            detect_complexity_structure(
                family_id=record.id,
                seq=str(record.seq),
            )
        )

    return calls

def local_complexity_metrics(seq: str, window: int = 40, step: int = 10) -> tuple[float, float, int]:
    seq = _clean_seq(seq)

    if len(seq) < window:
        return (
            shannon_entropy(seq),
            dominant_1mer_fraction(seq),
            max_homopolymer_run(seq),
        )

    min_entropy = 99.0
    max_dom1 = 0.0
    max_homo = 0

    for i in range(0, len(seq) - window + 1, step):
        w = seq[i:i + window]
        min_entropy = min(min_entropy, shannon_entropy(w))
        max_dom1 = max(max_dom1, dominant_1mer_fraction(w))
        max_homo = max(max_homo, max_homopolymer_run(w))

    return round(min_entropy, 3), round(max_dom1, 3), max_homo
