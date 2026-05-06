from .scoring import score_ltr, score_dna_tir, score_line
from .rules import is_ltr_candidate, is_dna_tir_candidate, is_line_candidate


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

    if not evidence:
        return "no_supporting_evidence"

    return ";".join(evidence)


def classify_family(f):
    candidates = []

    if is_ltr_candidate(f):
        candidates.append(("LTR", score_ltr(f)))

    if is_line_candidate(f):
        candidates.append(("LINE", score_line(f)))

    if is_dna_tir_candidate(f):
        candidates.append(("DNA_TIR", score_dna_tir(f)))

    if not candidates:
        return {
            "class": "Unknown_interspersed",
            "confidence": "LOW",
            "score": 0.0,
            "evidence": "no_class_rules_passed",
        }

    candidates.sort(key=lambda x: x[1], reverse=True)
    best_class, best_score = candidates[0]

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

    return {
        "class": best_class,
        "confidence": conf,
        "score": round(best_score, 3),
        "evidence": build_evidence_string(f),
    }
