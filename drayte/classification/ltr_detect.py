from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass
class LtrDetection:
    family_id: str
    consensus_len: int
    ltr_present: bool = False
    ltr_len: int = 0
    ltr_identity: float = 0.0
    left_ltr_start: int = 0
    left_ltr_end: int = 0
    right_ltr_start: int = 0
    right_ltr_end: int = 0
    internal_len: int = 0

    tg_ca_motif: bool = False
    left_ltr_motif: str = "NA"
    right_ltr_motif: str = "NA"
    pbs_like: bool = False
    ppt_like: bool = False

    ltr_score: float = 0.0
    ltr_structural_type: str = "none"
    confidence: str = "LOW"

def has_tg_ca_ltr_motif(seq: str, left_start: int, left_end: int, right_start: int, right_end: int) -> tuple[bool, str, str]:
    left = seq[left_start:left_end].upper()
    right = seq[right_start:right_end].upper()

    left_5 = left[:2]
    left_3 = left[-2:]
    right_5 = right[:2]
    right_3 = right[-2:]

    left_motif = f"{left_5}..{left_3}" if len(left) >= 4 else "NA"
    right_motif = f"{right_5}..{right_3}" if len(right) >= 4 else "NA"

    ok = (
        left_5 == "TG"
        and left_3 == "CA"
        and right_5 == "TG"
        and right_3 == "CA"
    )

    return ok, left_motif, right_motif

def has_pbs_like_region(seq: str, internal_start: int, window: int = 80) -> bool:
    region = seq[internal_start:internal_start + window].upper()
    if len(region) < 12:
        return False

    # Conservative PBS proxy: short, non-homopolymeric GC-balanced tract
    # near the 5' end of the internal region.
    for k in range(14, min(24, len(region)) + 1):
        for i in range(0, min(40, len(region) - k + 1)):
            sub = region[i:i + k]
            if "N" in sub:
                continue

            gc = (sub.count("G") + sub.count("C")) / k
            max_base = max(sub.count(b) for b in "ACGT") / k

            if 0.40 <= gc <= 0.65 and max_base <= 0.45:
                return True

    return False

def has_ppt_like_region(seq: str, internal_end: int, window: int = 120) -> bool:
    region = seq[max(0, internal_end - window):internal_end].upper()
    if len(region) < 10:
        return False

    # PPT: pyrimidine-rich tract upstream of 3' LTR.
    for k in range(10, min(30, len(region)) + 1):
        for i in range(0, len(region) - k + 1):
            sub = region[i:i + k]
            pyr = (sub.count("C") + sub.count("T")) / k
            if pyr >= 0.80 and sub.count("N") == 0:
                return True

    return False

def _identity(a: str, b: str) -> float:
    if not a or not b or len(a) != len(b):
        return 0.0

    matches = sum(x == y for x, y in zip(a, b))
    return matches / len(a)


def _clean_family_id(raw_id: str) -> str:
    return raw_id.split()[0].split("#")[0]

def infer_ltr_structural_type(
    consensus_len: int,
    ltr_len: int,
    internal_len: int,
) -> str:
    # Very small internal region: likely terminal/direct-repeat artifact,
    # not a credible LTR-retrotransposon-like structure.
    if internal_len < 100:
        return "none"

    ltr_fraction = ltr_len / consensus_len

    # reject implausibly tiny terminal repeats on large elements
    if consensus_len >= 5000 and ltr_fraction < 0.015:
        return "none"

    # TRIM
    if (
        consensus_len <= 1200
        and ltr_len <= 250
        and internal_len >= 100
        and ltr_fraction >= 0.02
    ):
        return "TRIM"

    # LARD-like
    if internal_len >= 1000 and consensus_len <= 7000:
        return "LARD_like"

    if consensus_len >= 1200 and ltr_len >= 80:
        return "LTR_like"

    return "none"

def detect_ltr_structure(
    family_id: str,
    seq: str,
    min_ltr_len: int = 80,
    max_ltr_len: int = 2000,
    min_identity: float = 0.75,
    terminal_slop: int = 30,
) -> LtrDetection:
    seq = (seq or "").upper()
    seq = "".join(x for x in seq if x in {"A", "C", "G", "T", "N"})
    consensus_len = len(seq)
    clean_id = _clean_family_id(family_id)

    result = LtrDetection(
        family_id=clean_id,
        consensus_len=consensus_len,
    )

    if consensus_len < (2 * min_ltr_len + 20):
        return result

    max_ltr_len = min(max_ltr_len, consensus_len // 2)
    best = None

    for ltr_len in range(min_ltr_len, max_ltr_len + 1):
        left_start = 0
        left_end = ltr_len
        right_start = consensus_len - ltr_len
        right_end = consensus_len

        left = seq[left_start:left_end]
        right = seq[right_start:right_end]

        ident = _identity(left, right)
        if ident < min_identity:
            continue

        if (
            best is None
            or ident > best["identity"]
            or (
                ident == best["identity"]
                and ltr_len > best["ltr_len"]
            )
        ):
            best = {
                "ltr_len": ltr_len,
                "identity": ident,
                "left_start": left_start,
                "left_end": left_end,
                "right_start": right_start,
                "right_end": right_end,
            }

    if best is None:
        return result

    ltr_len = best["ltr_len"]
    ltr_identity = best["identity"]
    left_start = best["left_start"]
    left_end = best["left_end"]
    right_start = best["right_start"]
    right_end = best["right_end"]
    internal_len = right_start - left_end

    ltr_structural_type = infer_ltr_structural_type(
        consensus_len=consensus_len,
        ltr_len=ltr_len,
        internal_len=internal_len,
    )

    if ltr_structural_type == "none":
        return result

    tg_ca_motif, left_ltr_motif, right_ltr_motif = has_tg_ca_ltr_motif(
        seq=seq,
        left_start=left_start,
        left_end=left_end,
        right_start=right_start,
        right_end=right_end,
    )

    pbs_like = has_pbs_like_region(
        seq=seq,
        internal_start=left_end,
    )

    ppt_like = has_ppt_like_region(
        seq=seq,
        internal_end=right_start,
    )

    ltr_score = score_ltr_structure(
        ltr_identity=ltr_identity,
        ltr_len=ltr_len,
        internal_len=internal_len,
        tg_ca_motif=tg_ca_motif,
        pbs_like=pbs_like,
        ppt_like=ppt_like,
    )

    if ltr_score >= 0.70:
        confidence = "HIGH"
    elif ltr_score >= 0.45:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"

    return LtrDetection(
        family_id=clean_id,
        consensus_len=consensus_len,
        ltr_present=True,
        ltr_len=ltr_len,
        ltr_identity=round(ltr_identity, 3),
        left_ltr_start=left_start,
        left_ltr_end=left_end,
        right_ltr_start=right_start,
        right_ltr_end=right_end,
        internal_len=internal_len,
        tg_ca_motif=tg_ca_motif,
        left_ltr_motif=left_ltr_motif,
        right_ltr_motif=right_ltr_motif,
        pbs_like=pbs_like,
        ppt_like=ppt_like,
        ltr_score=ltr_score,
        ltr_structural_type=ltr_structural_type,
        confidence=confidence,
    )

def score_ltr_structure(
    ltr_identity: float,
    ltr_len: int,
    internal_len: int,
    tg_ca_motif: bool = False,
    pbs_like: bool = False,
    ppt_like: bool = False,
) -> float:
    score = 0.0

    if ltr_identity >= 0.75:
        score += 0.25
    if ltr_identity >= 0.90:
        score += 0.20
    if ltr_len >= 80:
        score += 0.15
    if internal_len >= 100:
        score += 0.10
    if internal_len >= 1000:
        score += 0.10
    if tg_ca_motif:
        score += 0.10
# PBS proxy is currently diagnostic only; too permissive for scoring.
#    if pbs_like:
#        score += 0.05
    if ppt_like:
        score += 0.05

    return round(min(score, 1.0), 3)

def detect_ltrs_from_fasta(
    fasta: str | Path,
    min_ltr_len: int = 80,
    max_ltr_len: int = 2000,
    min_identity: float = 0.75,
) -> list[LtrDetection]:
    calls: list[LtrDetection] = []

    for record in SeqIO.parse(str(fasta), "fasta"):
        calls.append(
            detect_ltr_structure(
                family_id=record.id,
                seq=str(record.seq),
                min_ltr_len=min_ltr_len,
                max_ltr_len=max_ltr_len,
                min_identity=min_identity,
            )
        )

    return calls
