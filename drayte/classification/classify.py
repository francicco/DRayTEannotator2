from scoring import score_ltr, score_dna_tir, score_line
from rules import is_ltr_candidate, is_dna_tir_candidate, is_line_candidate


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
            "evidence": "no_class_rules_passed"
        }

    # pick best
    candidates.sort(key=lambda x: x[1], reverse=True)
    best_class, best_score = candidates[0]

    # margin
    if len(candidates) > 1:
        margin = best_score - candidates[1][1]
    else:
        margin = best_score

    # confidence
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
        "evidence": build_evidence_string(f)
    }
