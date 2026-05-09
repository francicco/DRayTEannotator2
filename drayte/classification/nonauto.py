from __future__ import annotations

from dataclasses import dataclass

from .models import Family


@dataclass
class NonAutonomousCall:
    family_id: str
    candidate_type: str
    score: float
    confidence: str
    evidence: str


def has_coding_domains(f: Family) -> bool:
    return bool(
        f.rt_present
        or f.integrase_present
        or f.transposase_present
        or f.rnaseh_present
        or f.gag_present
        or f.helitron_domain_present
    )


def has_substantial_orf(
    f: Family,
    max_orf_nt: int = 300,
) -> bool:
    return f.orf_max_len >= max_orf_nt


def is_short_element(
    f: Family,
    max_len: int,
) -> bool:
    return f.consensus_len > 0 and f.consensus_len <= max_len


def score_sine_candidate(f: Family) -> NonAutonomousCall | None:
    """
    Conservative SINE candidate rule.

    This does not prove SINE status. It identifies short, non-coding,
    polyA-associated retroelement-like families for later validation.
    """

    evidence = []
    score = 0.0

    if not is_short_element(f, 700):
        return None

    evidence.append("short_len<=700")
    score += 0.25

    if has_substantial_orf(f):
        return None

    evidence.append("no_substantial_orf")
    score += 0.25

    if has_coding_domains(f):
        return None

    evidence.append("no_coding_domains")
    score += 0.20

    if f.polyA_present:
        evidence.append("polyA_present")
        score += 0.25
    else:
        return None

    if f.homology_order == "SINE" or f.dfam_order == "SINE":
        evidence.append("sine_homology")
        score += 0.30

    if f.rescue_order == "SINE":
        evidence.append("sine_rescue")
        score += 0.20

    score = min(score, 0.99)

    confidence = "LOW"
    if score >= 0.75:
        confidence = "MEDIUM"
    if score >= 0.90:
        confidence = "HIGH"

    return NonAutonomousCall(
        family_id=f.family_id,
        candidate_type="SINE_candidate",
        score=round(score, 3),
        confidence=confidence,
        evidence=";".join(evidence),
    )


def score_mite_candidate(f: Family) -> NonAutonomousCall | None:
    """
    Conservative MITE candidate rule.

    Baseline MITE evidence is short + TIR + non-coding.
    HIGH confidence requires TSD support or TIR homology/rescue.
    """

    evidence = []
    score = 0.0

    if not is_short_element(f, 800):
        return None

    evidence.append("short_len<=800")
    score += 0.15

    if not f.tir_present:
        return None

    evidence.append("tir_present")
    score += 0.30

    if has_substantial_orf(f):
        return None

    evidence.append("no_substantial_orf")
    score += 0.20

    if f.transposase_present:
        return None

    evidence.append("no_transposase_domain")
    score += 0.10

    if f.tsd_present:
        evidence.append("tsd_present")
        score += 0.20

    if (
        f.homology_order == "TIR"
        or f.dfam_order == "TIR"
        or f.rescue_order == "TIR"
    ):
        evidence.append("tir_partner_or_homology")
        score += 0.20

    score = min(score, 0.99)

    confidence = "LOW"

    if score >= 0.70:
        confidence = "MEDIUM"

    if score >= 0.90:
        confidence = "HIGH"

    return NonAutonomousCall(
        family_id=f.family_id,
        candidate_type="MITE_candidate",
        score=round(score, 3),
        confidence=confidence,
        evidence=";".join(evidence),
    )

def score_trim_candidate(f: Family) -> NonAutonomousCall | None:
    """
    TRIM candidate: short LTR-like non-autonomous element.
    """

    evidence = []
    score = 0.0

    if not is_short_element(f, 1500):
        return None

    evidence.append("short_len<=1500")
    score += 0.20

    if not f.ltr_present:
        return None

    evidence.append("ltr_present")
    score += 0.35

    if f.rt_present or f.integrase_present:
        return None

    evidence.append("no_rt_integrase")
    score += 0.25

    if not has_substantial_orf(f):
        evidence.append("no_substantial_orf")
        score += 0.20

    score = min(score, 0.99)

    confidence = "LOW"
    if score >= 0.70:
        confidence = "MEDIUM"
    if score >= 0.90:
        confidence = "HIGH"

    return NonAutonomousCall(
        family_id=f.family_id,
        candidate_type="TRIM_candidate",
        score=round(score, 3),
        confidence=confidence,
        evidence=";".join(evidence),
    )


def infer_nonautonomous_candidates(
    families: list[Family],
) -> list[NonAutonomousCall]:
    calls: list[NonAutonomousCall] = []

    for f in families:
        for detector in (
            score_sine_candidate,
            score_mite_candidate,
            score_trim_candidate,
        ):
            call = detector(f)

            if call is not None:
                calls.append(call)

    return calls
