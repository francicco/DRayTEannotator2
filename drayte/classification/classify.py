
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

rescue_label_map = {
    "DNA": "DNA_TIR",
    "TIR": "DNA_TIR",
    "Helitron": "HELITRON",
    "LINE": "LINE",
    "LTR": "LTR",
    "SINE": "SINE",
}

def build_evidence_string(f):
    evidence = []

    if f.ltr_present:
        evidence.append(f"LTR_structure={f.ltr_structural_type}")

    if f.tg_ca_motif:
        evidence.append("TG_CA_LTR_motif")

    if f.ppt_like:
        evidence.append("PPT_like")

    if (
        f.ltr_present
        and (
            f.homology_order == "LINE"
            or f.dfam_order == "LINE"
            or f.header_class == "LINE"
            or f.header_superfamily == "Penelope"
            or f.homology_superfamily == "Penelope"
            or f.dfam_superfamily == "Penelope"
        )
    ):
        evidence.append(
            "CONFLICT:LTR_structure_with_LINE_or_Penelope_signal"
        )

    if f.rnaseh_present:
        evidence.append("RNaseH_domain")

    if f.gag_present:
        evidence.append("GAG_domain")

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
        evidence.append(
            f"TSD:{f.tsd_len}:{f.tsd_seq}:support={f.tsd_support}"
        )

    if f.structural_superfamily_hint not in {
        "", "NA", "Unknown", "unknown", None
    }:
        evidence.append(
            f"structural_superfamily={f.structural_superfamily_hint}:"
            f"confidence={f.structural_superfamily_confidence}"
        )

    if f.polyA_present:
        evidence.append("polyA")

    if f.dfam_model:
        evidence.append(
            f"dfam={f.dfam_model}:{f.dfam_order}:"
            f"{f.dfam_superfamily}:{f.dfam_score}"
        )

    if f.rescue_superfamily not in {"", "NA", "Unknown", None}:
        evidence.append(
            f"rescue={f.rescue_target}:{f.rescue_order}:"
            f"{f.rescue_superfamily}:"
            f"id={f.rescue_identity}:"
            f"len={f.rescue_aln_len}:"
            f"bits={f.rescue_bits}"
        )

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

def collapse_candidates(candidates):
    best = {}

    for label, score in candidates:
        if label not in best or score > best[label]:
            best[label] = score

    return list(best.items())

def compatible_superfamily(order: str, superfamily: str) -> bool:
    if superfamily in {"", "NA", "Unknown", None}:
        return False

    sf = str(superfamily)

    if order == "Helitron":
        return "Helitron" in sf

    if order == "TIR":
        sf_lower = sf.lower()
        return (
            "helitron" not in sf_lower
            and sf not in {
                "Gypsy", "Copia", "Bel-Pao", "Pao", "DIRS", "ERV1",
                "Penelope", "L1", "L2", "R1", "R2", "RTE", "I",
                "Jockey", "CRE",
            }
        )

    if order == "LTR":
        return (
            sf in {"Gypsy", "Copia", "Bel-Pao", "Pao", "DIRS", "ERV1"}
            or "LTR" in sf
        )

    if order == "LINE":
        return (
            "Helitron" not in sf
            and sf not in {"Gypsy", "Copia", "Bel-Pao", "Pao", "DIRS", "ERV1"}
        )

    if order == "SINE":
        return "Helitron" not in sf

    return True

def count_supporting_evidence(f, best_label):
    support = 0

    if f.dfam_order not in {"", "Unknown", None}:
        support += 1

    if best_label == "LTR" and f.ltr_present:
        support += 1

    if best_label == "DNA_TIR" and f.tir_present:
        support += 1

    if best_label == "HELITRON" and f.helitron_signal:
        support += 1

    if best_label in {"LINE", "LTR"} and f.rt_present:
        support += 1

    if best_label == "LTR" and f.integrase_present:
        support += 1

    if best_label == "DNA_TIR" and f.transposase_present:
        support += 1

    if f.rescue_order not in {"", "Unknown", None}:
        support += 1

    if f.homology_order not in {"", "Unknown", None} and f.homology_score >= 0.8:
        support += 1

    return support

def classify_family(f):
    candidates = []

    homology_order_label_map = {
        "DNA": "DNA_TIR",
        "TIR": "DNA_TIR",
        "Helitron": "HELITRON",
        "SINE": "SINE",
        "LTR": "LTR",
        "LINE": "LINE",
    }

    if (
        f.homology_score >= 0.8
        and f.homology_order in homology_order_label_map
    ):
        candidates.append(
            (
                homology_order_label_map[f.homology_order],
                min(0.99, max(0.85, f.homology_score)),
            )
        )

    # Hard prior: trusted HELIANO / header annotation
    if (
        f.header_class in {"RC", "Helitron"}
        or (
            f.header_superfamily not in {"", "Unknown", "unknown", "NA", None}
            and "helitron" in f.header_superfamily.lower()
        )
    ):
        candidates.append(("HELITRON", 0.95))

    if f.dfam_order == "LTR":
        candidates.append(("LTR", min(0.99, max(0.85, f.dfam_score / 100.0))))

    if f.dfam_order == "LINE":
        candidates.append(("LINE", min(0.99, max(0.85, f.dfam_score / 100.0))))

    if f.dfam_order == "TIR":
        candidates.append(("DNA_TIR", min(0.99, max(0.85, f.dfam_score / 100.0))))

    if f.dfam_order == "Helitron":
        candidates.append(("HELITRON", min(0.99, max(0.90, f.dfam_score / 100.0))))

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

    rescue_candidates = []

    rescue_label = rescue_label_map.get(
        f.rescue_order
    )

    if (
        rescue_label is not None
        and f.rescue_identity >= 0.80
    ):
        rescue_candidates.append(
            (
                rescue_label,
                min(0.75, f.rescue_identity),
            )
        )

    if not candidates and rescue_candidates:
        candidates = rescue_candidates

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
    if f.homology_score >= 0.8 and f.homology_order in homology_order_label_map:
        preferred = homology_order_label_map[f.homology_order]
        candidates = [
            c for c in candidates
            if c[0] == preferred
        ] or candidates

    candidates = collapse_candidates(candidates)

    candidates.sort(key=lambda x: x[1], reverse=True)
    best_label, best_score = candidates[0]

    if len(candidates) > 1:
        margin = best_score - candidates[1][1]
    else:
        margin = best_score

    support_count = count_supporting_evidence(f, best_label)

    if (
        best_score >= 0.80
        and margin >= 0.15
        and support_count >= 2
    ):
        conf = "HIGH"
    elif best_score >= 0.60 and margin >= 0.10:
        conf = "MEDIUM"
    else:
        conf = "LOW"        

    te_class, order = CLASS_MAP[best_label]

    superfamily = "Unknown"

    for candidate_superfamily in [
        f.homology_superfamily,
        f.dfam_superfamily,
        f.rescue_superfamily,
    ]:
        if compatible_superfamily(order, candidate_superfamily):
            superfamily = candidate_superfamily
            break

    if (
        superfamily == "Unknown"
        and best_label == "DNA_TIR"
        and f.structural_superfamily_hint not in {"", "NA", "Unknown", "unknown", None}
        and compatible_superfamily(order, f.structural_superfamily_hint)
    ):
        superfamily = f.structural_superfamily_hint

    return {
        "class": te_class,
        "order": order,
        "superfamily": superfamily,
        "status": infer_status(candidates, best_score, margin),
        "confidence": conf,
        "score": round(best_score, 3),
        "evidence": build_evidence_string(f),
    }
