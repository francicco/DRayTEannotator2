from .scoring import score_ltr, score_dna_tir, score_line
from .rules import is_ltr_candidate, is_dna_tir_candidate, is_line_candidate

from .scoring import (
    score_ltr,
    score_dna_tir,
    score_line,
    score_helitron,
    score_sine,
)

from .rules import (
    is_ltr_candidate,
    is_dna_tir_candidate,
    is_line_candidate,
    is_helitron_candidate,
    is_sine_candidate,
)

CLASS_MAP = {
    "LTR": ("Class_I", "LTR"),
    "LINE": ("Class_I", "LINE"),
    "DNA_TIR": ("Class_II", "TIR"),
    "HELITRON": ("Class_II", "Helitron"),
    "SINE": ("Class_I", "SINE"),
}


def build_evidence_string(f):
    evidence = []

    if f.ltr_present:
        evidence.append("LTR_structure")
    if f.tir_present:
        evidence.append("TIR_structure")
    if f.helitron_signal:
        evidence.append("Helitron_signal")
    if f.rt_present:
        evidence.append("RT_domain")
    if f.integrase_present:
        evidence.append("Integrase_domain")
    if f.transposase_present:
        evidence.append("Transposase_domain")
    if f.tsd_present:
        evidence.append("TSD")
    if f.polyA_present:
        evidence.append("polyA")
    if f.homology_class not in {"", "NA", "Unknown", None}:
        evidence.append(f"homology={f.homology_class}:{f.homology_superfamily}:{f.homology_score}")

    return ";".join(evidence) if evidence else "no_supporting_evidence"


def infer_status(candidates, best_score, margin):
    if not candidates:
        return "unknown"

    if len(candidates) > 1 and margin < 0.10:
        return "ambiguous"

    if len(candidates) > 1:
        return "conflicting_evidence"

    if best_score < 0.50:
        return "weak_evidence"

    return "OK"


def classify_family(f):
    candidates = []

    if f.dfam_order == "LTR":
        candidates.append(("LTR", max(0.85, f.dfam_score / 50.0)))
    
    if f.dfam_order == "LINE":
        candidates.append(("LINE", max(0.85, f.dfam_score / 50.0)))
    
    if f.dfam_order == "TIR":
        candidates.append(("DNA_TIR", max(0.85, f.dfam_score / 50.0)))
    
    if f.dfam_order == "Helitron":
        candidates.append(("HELITRON", max(0.90, f.dfam_score / 50.0)))

    if is_ltr_candidate(f):
        candidates.append(("LTR", score_ltr(f)))

    if is_line_candidate(f):
        candidates.append(("LINE", score_line(f)))

    if is_dna_tir_candidate(f):
        candidates.append(("DNA_TIR", score_dna_tir(f)))

    if is_helitron_candidate(f):
        candidates.append(("HELITRON", score_helitron(f)))

    if is_sine_candidate(f):
        candidates.append(("SINE", score_sine(f)))

    if not candidates:
        return {
            "class": "Unknown",
            "order": "Unknown",
            "superfamily": "Unknown",
            "status": "unknown",
            "confidence": "LOW",
            "score": 0.0,
            "evidence": "no_class_rules_passed",
        }

    # High-confidence RepeatMasker-style homology priors should not be
    # overridden by generic protein domains such as DDE or RT-like hits.
    homology_label_map = {
        "DNA": "DNA_TIR",
        "Helitron": "HELITRON",
        "SINE": "SINE",
        "LTR": "LTR",
        "LINE": "LINE",
    }

    if f.homology_score >= 0.8 and f.homology_class in homology_label_map:
        preferred = homology_label_map[f.homology_class]
        candidates = [
            c for c in candidates
            if c[0] == preferred
        ] or candidates

    candidates.sort(key=lambda x: x[1], reverse=True)
    best_label, best_score = candidates[0]

    if len(candidates) > 1:
        margin = best_score - candidates[1][1]
    else:
        margin = best_score

    if best_score >= 0.80 and margin >= 0.15:
        conf = "HIGH"
    elif best_score >= 0.60 and margin >= 0.10:
        conf = "MEDIUM"
    else:
        conf = "LOW"

    te_class, order = CLASS_MAP[best_label]

    if f.dfam_superfamily not in {"", "NA", "Unknown", None}:
        superfamily = f.dfam_superfamily
    else:
        superfamily = f.homology_superfamily

    if f.dfam_model:
        evidence.append(
            f"dfam={f.dfam_model}:{f.dfam_order}:{f.dfam_superfamily}:{f.dfam_score}"
        )

    if superfamily in {"", "NA", None}:
        superfamily = "Unknown"

    return {
        "class": te_class,
        "order": order,
        "superfamily": superfamily,
        "status": infer_status(candidates, best_score, margin),
        "confidence": conf,
        "score": round(best_score, 3),
        "evidence": build_evidence_string(f),
    }
