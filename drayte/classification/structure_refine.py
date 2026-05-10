from __future__ import annotations

from dataclasses import dataclass


LOW_COMPLEXITY_KMERS = {
    "AT",
    "TA",
    "CG",
    "GC",
    "AA",
    "TT",
    "CC",
    "GG",
}


@dataclass
class RefinedStructureCall:
    family_id: str
    structure_class: str
    confidence: str
    refined_score: float
    superfamily_hint: str = "NA"
    notes: str = ""


# ---------------------------------------------------------
# Complexity
# ---------------------------------------------------------


def shannon_complexity(seq: str) -> float:
    seq = seq.upper()

    if not seq:
        return 0.0

    freqs = {}

    for b in seq:
        freqs[b] = freqs.get(b, 0) + 1

    import math

    h = 0.0

    for n in freqs.values():
        p = n / len(seq)
        h -= p * math.log2(p)

    return h



def tir_is_low_complexity(seq: str) -> bool:
    seq = seq.upper()

    if len(set(seq)) <= 2:
        return True

    if seq in LOW_COMPLEXITY_KMERS:
        return True

    if shannon_complexity(seq) < 1.2:
        return True

    return False


# ---------------------------------------------------------
# TSD-superfamily heuristics
# ---------------------------------------------------------


def infer_tir_superfamily(tsd_seq: str, tsd_len: int) -> str:
    tsd_seq = (tsd_seq or "").upper()

    if tsd_seq == "TTAA":
        return "piggyBac"

    if tsd_seq in {"TA", "TAA", "TTA"}:
        return "Tc1/Mariner"

    if 8 <= tsd_len <= 11:
        return "PIF/Harbinger"

    if tsd_len == 9:
        return "hAT"

    return "NA"


# ---------------------------------------------------------
# TIR quality refinement
# ---------------------------------------------------------


def refine_tir_detection(det):
    if not det.tir_present:
        return det

    left = det.left_tir_seq.upper()
    right = det.right_tir_seq.upper()

    if tir_is_low_complexity(left):
        det.tir_present = False
        det.confidence = "LOW"
        return det

    if tir_is_low_complexity(right):
        det.tir_present = False
        det.confidence = "LOW"
        return det

    if det.tir_len < 10:
        det.confidence = "LOW"

    return det


# ---------------------------------------------------------
# Arbitration logic
# ---------------------------------------------------------


def arbitrate_nonautonomous_calls(calls):
    by_family = {}

    for c in calls:
        fam = c.family_id

        if fam not in by_family:
            by_family[fam] = []

        by_family[fam].append(c)

    final = []

    for fam, fam_calls in by_family.items():
        if len(fam_calls) == 1:
            final.extend(fam_calls)
            continue

        mite = [x for x in fam_calls if x.structure_class == "MITE"]
        sine = [x for x in fam_calls if x.structure_class == "SINE"]

        if mite and sine:
            best_mite = sorted(
                mite,
                key=lambda x: x.structural_score,
                reverse=True,
            )[0]

            final.append(best_mite)
            continue

        best = sorted(
            fam_calls,
            key=lambda x: x.structural_score,
            reverse=True,
        )[0]

        final.append(best)

    return final


# ---------------------------------------------------------
# MITE refinement
# ---------------------------------------------------------


def refine_mite_call(call):
    if call.structure_class != "MITE":
        return call

    hint = infer_tir_superfamily(
        call.tsd_seq,
        call.tsd_len,
    )

    call.superfamily_hint = hint

    if hint != "NA":
        call.structural_score += 0.15

    return call
