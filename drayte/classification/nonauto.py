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


def confidence_from_score(score: float) -> str:
    if score >= 0.90:
        return "HIGH"

    if score >= 0.70:
        return "MEDIUM"

    return "LOW"


def tir_grade(f: Family) -> str:
    """
    Convert raw TIR evidence into a biological confidence grade.

    This is deliberately stricter than a binary TIR flag. Short low-complexity
    inverted repeats can occur by chance in short consensuses, so MITE calling
    should prefer terminal, sufficiently long, high-identity TIRs.
    """

    if not f.tir_present:
        return "ABSENT"

    conf = (getattr(f, "tir_confidence", "LOW") or "LOW").upper()
    tir_len = int(getattr(f, "tir_len", 0) or 0)
    tir_identity = float(getattr(f, "tir_identity", 0.0) or 0.0)

    if conf == "HIGH" and tir_len >= 16 and tir_identity >= 0.90:
        return "STRONG"

    if conf in {"HIGH", "MEDIUM"} and tir_len >= 14 and tir_identity >= 0.85:
        return "MODERATE"

    return "WEAK"


def infer_mite_superfamily_from_tsd(f: Family) -> str:
    """
    Infer only a structural superfamily hint from TSD pattern.

    This should not be treated as a definitive taxonomic assignment unless
    supported by homology/domain evidence. TSD length is useful for TIR
    elements, but several superfamilies overlap.
    """

    tsd_seq = (getattr(f, "tsd_seq", "") or "").upper()
    tsd_len = int(getattr(f, "tsd_len", 0) or 0)

    if tsd_seq == "TA":
        return "TcMar-Mariner"

    if tsd_seq == "TTAA":
        return "piggyBac_like"

    if tsd_seq == "TAA" or tsd_len == 3:
        return "PIF-Harbinger_or_CACTA_like"

    if tsd_len == 8:
        return "hAT_or_P_or_Merlin_like"

    if tsd_len in {9, 10, 11}:
        return "Mutator_like"

    if tsd_len in {2, 3}:
        return "CACTA_or_PIF-Harbinger_like"

    return "Unknown"


def score_sine_candidate(f: Family) -> NonAutonomousCall | None:
    """
    Conservative SINE candidate rule.

    SINEs are short non-autonomous class-I elements derived from Pol-III RNAs.
    They should not have terminal inverted repeats. A TIR call is therefore a
    hard exclusion, not merely competing evidence.
    """

    evidence: list[str] = []
    score = 0.0

    if f.tir_present:
        return None

    if f.ltr_present:
        return None

    if not is_short_element(f, 700):
        return None

    evidence.append("short_len<=700")
    score += 0.20

    if has_substantial_orf(f):
        return None

    evidence.append("no_substantial_orf")
    score += 0.20

    if has_coding_domains(f):
        return None

    evidence.append("no_coding_domains")
    score += 0.20

    if not f.polyA_present:
        return None

    evidence.append("polyA_present")
    score += 0.25

    if f.tsd_present:
        # SINEs can have TSDs, but the length/sequence is partner dependent and
        # is not enough to classify a SINE by itself.
        evidence.append(f"tsd_len={int(getattr(f, 'tsd_len', 0) or 0)}")
        score += 0.05

    if f.homology_order == "SINE" or f.dfam_order == "SINE":
        evidence.append("sine_homology")
        score += 0.30

    if f.rescue_order == "SINE":
        evidence.append("sine_rescue")
        score += 0.20

    score = min(score, 0.99)

    return NonAutonomousCall(
        family_id=f.family_id,
        candidate_type="SINE_candidate",
        score=round(score, 3),
        confidence=confidence_from_score(score),
        evidence=";".join(evidence),
    )


def score_mite_candidate(f: Family) -> NonAutonomousCall | None:
    """
    Conservative MITE candidate rule.

    A MITE is a short, non-autonomous TIR element. TSD evidence is not required
    for every family because consensus boundaries may be imperfect, but a weak
    TIR without TSD support should not pass.
    """

    evidence: list[str] = []
    score = 0.0

    if not is_short_element(f, 800):
        return None

    evidence.append("short_len<=800")
    score += 0.15

    grade = tir_grade(f)
    if grade == "ABSENT":
        return None

    if grade == "WEAK" and not f.tsd_present:
        return None

    evidence.append(f"tir_grade={grade}")

    if grade == "STRONG":
        score += 0.35
    elif grade == "MODERATE":
        score += 0.25
    else:
        score += 0.10

    if has_substantial_orf(f):
        return None

    evidence.append("no_substantial_orf")
    score += 0.20

    if f.transposase_present:
        return None

    evidence.append("no_transposase_domain")
    score += 0.10

    if f.rt_present or f.integrase_present or f.rnaseh_present or f.gag_present:
        return None

    evidence.append("no_retroelement_domains")
    score += 0.10

    if f.tsd_present:
        tsd_len = int(getattr(f, "tsd_len", 0) or 0)
        tsd_seq = (getattr(f, "tsd_seq", "") or "NA").upper()
        tsd_support = float(getattr(f, "tsd_support", 0.0) or 0.0)

        evidence.append("tsd_present")
        evidence.append(f"tsd_len={tsd_len}")
        evidence.append(f"tsd_seq={tsd_seq}")
        evidence.append(f"tsd_support={tsd_support}")

        inferred_sf = infer_mite_superfamily_from_tsd(f)

        if inferred_sf != "Unknown":
            evidence.append(f"tsd_superfamily_hint={inferred_sf}")
            score += 0.20
        else:
            score += 0.10

        if tsd_support >= 0.30:
            evidence.append("tsd_support>=0.30")
            score += 0.05

        if tsd_support >= 0.50:
            evidence.append("tsd_support>=0.50")
            score += 0.05

    if (
        f.homology_order == "TIR"
        or f.dfam_order == "TIR"
        or f.rescue_order == "TIR"
    ):
        evidence.append("tir_partner_or_homology")
        score += 0.20

    score = min(score, 0.99)

    return NonAutonomousCall(
        family_id=f.family_id,
        candidate_type="MITE_candidate",
        score=round(score, 3),
        confidence=confidence_from_score(score),
        evidence=";".join(evidence),
    )


def _score_ltr_like_nonautonomous(
    f: Family,
    candidate_type: str,
    min_len: int,
    max_len: int | None,
) -> NonAutonomousCall | None:
    evidence: list[str] = []
    score = 0.0

    if f.consensus_len <= 0:
        return None

    if f.consensus_len < min_len:
        return None

    if max_len is not None and f.consensus_len > max_len:
        return None

    if max_len is None:
        evidence.append(f"len>={min_len}")
    else:
        evidence.append(f"{min_len}<=len<={max_len}")
    score += 0.15

    if not f.ltr_present:
        return None

    evidence.append("ltr_present")
    score += 0.35

    # Non-autonomous LTR derivatives should lack the enzymatic retroelement core.
    if f.rt_present or f.integrase_present or f.rnaseh_present:
        return None

    evidence.append("no_rt_integrase_rnaseh")
    score += 0.25

    if has_substantial_orf(f):
        return None

    evidence.append("no_substantial_orf")
    score += 0.20

    if f.gag_present:
        evidence.append("gag_only_or_partial_gag")
        score += 0.05

    score = min(score, 0.99)

    return NonAutonomousCall(
        family_id=f.family_id,
        candidate_type=candidate_type,
        score=round(score, 3),
        confidence=confidence_from_score(score),
        evidence=";".join(evidence),
    )


def score_trim_candidate(f: Family) -> NonAutonomousCall | None:
    """TRIM candidate: short LTR-like non-autonomous element."""

    return _score_ltr_like_nonautonomous(
        f=f,
        candidate_type="TRIM_candidate",
        min_len=80,
        max_len=1500,
    )


def score_lard_candidate(f: Family) -> NonAutonomousCall | None:
    """LARD candidate: longer LTR-like non-autonomous element."""

    return _score_ltr_like_nonautonomous(
        f=f,
        candidate_type="LARD_candidate",
        min_len=1501,
        max_len=None,
    )


def arbitrate_nonautonomous_calls(
    calls: list[NonAutonomousCall],
) -> list[NonAutonomousCall]:
    """
    Keep one non-autonomous call per family.

    Higher score wins. On ties, prefer structurally diagnostic classes.
    """

    priority = {
        "MITE_candidate": 4,
        "TRIM_candidate": 3,
        "LARD_candidate": 3,
        "SINE_candidate": 2,
    }

    best_by_family: dict[str, NonAutonomousCall] = {}

    for call in calls:
        current = best_by_family.get(call.family_id)

        if current is None:
            best_by_family[call.family_id] = call
            continue

        if call.score > current.score:
            best_by_family[call.family_id] = call
            continue

        if call.score == current.score:
            if priority.get(call.candidate_type, 0) > priority.get(
                current.candidate_type,
                0,
            ):
                best_by_family[call.family_id] = call

    return [
        best_by_family[k]
        for k in sorted(best_by_family)
    ]


def infer_nonautonomous_candidates(
    families: list[Family],
) -> list[NonAutonomousCall]:
    calls: list[NonAutonomousCall] = []

    for f in families:
        for detector in (
            score_sine_candidate,
            score_mite_candidate,
            score_trim_candidate,
            score_lard_candidate,
        ):
            call = detector(f)

            if call is not None:
                calls.append(call)

    return arbitrate_nonautonomous_calls(calls)

