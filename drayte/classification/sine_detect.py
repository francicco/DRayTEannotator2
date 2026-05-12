from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from .ids import clean_family_id


@dataclass
class SineDetection:
    family_id: str
    consensus_len: int = 0

    polyA_present: bool = False
    polyA_tail_len: int = 0
    polyA_fraction: float = 0.0
    tail_type: str = "none"

    poliii_a_box: bool = False
    poliii_b_box: bool = False
    poliii_score: float = 0.0

    sine_head_type: str = "unknown"
    sine_head_score: float = 0.0

    sine_score: float = 0.0
    sine_candidate: bool = False
    confidence: str = "LOW"


def _longest_run(seq: str, base: str) -> int:
    best = 0
    cur = 0
    for c in seq.upper():
        if c == base:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def _at_fraction(seq: str) -> float:
    cleaned = "".join(c for c in seq.upper() if c in {"A", "C", "G", "T"})
    if not cleaned:
        return 0.0
    return (cleaned.count("A") + cleaned.count("T")) / len(cleaned)


def _terminal_run(seq: str, base: str) -> int:
    n = 0
    for c in reversed(seq.upper()):
        if c == base:
            n += 1
        else:
            break
    return n


# Deliberately broad Pol III A/B-box motifs. They are weak evidence only;
# confirmed SINE subclassification should later use tRNA/5S/7SL covariance models.
_A_BOX_RE = re.compile(r"[TC][AG]G[CT][ACGT]{2}A[AG][ACGT]G")
_B_BOX_RE = re.compile(r"GGTT[CG]GA[ACGT]{1,3}CC")


def _infer_tail_type(seq: str, poly_run: int, at_fraction: float) -> str:
    seq = (seq or "").upper()

    longest_a = _longest_run(seq, "A")
    longest_t = _longest_run(seq, "T")
    terminal_a = _terminal_run(seq, "A")
    terminal_t = _terminal_run(seq, "T")

    if max(longest_a, terminal_a) >= 10 and max(longest_a, terminal_a) >= max(longest_t, terminal_t):
        return "polyA"

    if max(longest_t, terminal_t) >= 10:
        return "polyT"

    if at_fraction >= 0.70:
        a_count = seq.count("A")
        t_count = seq.count("T")

        if a_count > t_count * 1.5:
            return "A_rich"

        if t_count > a_count * 1.5:
            return "T_rich"

        return "mixed_AT"

    return "none"

def detect_poliii_boxes(head: str) -> tuple[bool, bool, float]:
    head = (head or "").upper()
    a_box = bool(_A_BOX_RE.search(head))
    b_box = bool(_B_BOX_RE.search(head))

    score = 0.0
    if a_box:
        score += 0.5
    if b_box:
        score += 0.5

    return a_box, b_box, score


def detect_sine_structure(
    family_id: str,
    seq: str,
    max_len: int = 700,
    tail_window: int = 80,
    head_window: int = 120,
    min_poly_run: int = 10,
    min_tail_at_fraction: float = 0.65,
    min_score: float = 0.75,
) -> SineDetection:
    seq = (seq or "").upper()
    seq = "".join(c for c in seq if c in {"A", "C", "G", "T", "N"})
    consensus_len = len(seq)

    tail = seq[-tail_window:] if seq else ""
    head = seq[:head_window] if seq else ""

    longest_a = _longest_run(tail, "A")
    longest_t = _longest_run(tail, "T")
    terminal_a = _terminal_run(tail, "A")
    terminal_t = _terminal_run(tail, "T")
    poly_run = max(longest_a, longest_t, terminal_a, terminal_t)

    # Use the strongest A/T enrichment in the broader 80-bp tail or the
    # immediate terminal 40 bp. Consensus tails are often short enough that
    # the broader window contains unrelated internal sequence.
    terminal_tail = tail[-40:]
    polyA_fraction = max(_at_fraction(tail), _at_fraction(terminal_tail))
    tail_type = _infer_tail_type(
        terminal_tail,
        poly_run=poly_run,
        at_fraction=polyA_fraction,
    )
    polyA_present = bool(
        poly_run >= min_poly_run
        and tail_type in {"polyA", "polyT"}
        and (
            polyA_fraction >= min_tail_at_fraction
            or max(terminal_a, terminal_t) >= min_poly_run
        )
    )

    poliii_a_box, poliii_b_box, poliii_score = detect_poliii_boxes(head)
    sine_head_type, sine_head_score = detect_sine_head(head)

    score = 0.0
    if consensus_len and consensus_len <= max_len:
        score += 0.30
    if polyA_present:
        score += 0.40
    elif tail_type in {"A_rich", "T_rich", "mixed_AT"}:
        score += 0.15
    
    if polyA_fraction >= 0.70:
        score += 0.15    
    if poliii_score > 0:
        score += 0.10

    score = round(min(score, 1.0), 3)
    sine_candidate = bool(
        score >= min_score
        and (
            poly_run >= 8
            or tail_type in {"polyA", "polyT"}
        )
    )

    if sine_head_score >= 0.7:
        score += 0.20
    elif sine_head_score > 0:
        score += 0.10

    if score >= 0.85:
        confidence = "HIGH"
    elif score >= 0.70:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"

    return SineDetection(
        family_id=clean_family_id(family_id),
        consensus_len=consensus_len,
        polyA_present=polyA_present,
        polyA_tail_len=poly_run,
        polyA_fraction=round(polyA_fraction, 3),
        tail_type=tail_type,
        poliii_a_box=poliii_a_box,
        poliii_b_box=poliii_b_box,
        poliii_score=round(poliii_score, 3),
        sine_head_type=sine_head_type,
        sine_head_score=round(sine_head_score, 3),
        sine_score=score,
        sine_candidate=sine_candidate,
        confidence=confidence,
    )

def detect_sine_head(head: str) -> tuple[str, float]:
    head = (head or "").upper()

    # Very lightweight heuristics.
    # These are not replacements for Infernal/Rfam; they are confidence hints.

    a_box_like = bool(re.search(r"[AT]GG[CT][ATGC]{2,6}A[ATGC]{2,6}G", head))
    b_box_like = bool(re.search(r"GGTTC[AG][ATGC]{2,8}CC", head))

    if a_box_like and b_box_like:
        return "tRNA_like", 0.7

    if a_box_like or b_box_like:
        return "PolIII_like", 0.4

    # 5S-derived SINEs can be GC-rich in the head, but this is weak.
    head_80 = head[:80]
    if len(head_80) >= 40:
        gc = (head_80.count("G") + head_80.count("C")) / len(head_80)
        if gc >= 0.60:
            return "possible_5S_like", 0.25

    return "unknown", 0.0

def detect_sines_from_fasta(
    fasta: str | Path,
    max_len: int = 700,
    tail_window: int = 80,
    head_window: int = 120,
    min_poly_run: int = 10,
    min_tail_at_fraction: float = 0.65,
    min_score: float = 0.75,
) -> list[SineDetection]:
    calls: list[SineDetection] = []

    for rec in SeqIO.parse(str(fasta), "fasta"):
        calls.append(
            detect_sine_structure(
                family_id=rec.id,
                seq=str(rec.seq),
                max_len=max_len,
                tail_window=tail_window,
                head_window=head_window,
                min_poly_run=min_poly_run,
                min_tail_at_fraction=min_tail_at_fraction,
                min_score=min_score,
            )
        )

    return calls
