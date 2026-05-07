from __future__ import annotations

import re


DFAM_MODEL_RULES = [
    (r"helitron", ("Class_II", "Helitron", "Helitron")),
    (r"gypsy", ("Class_I", "LTR", "Gypsy")),
    (r"copia", ("Class_I", "LTR", "Copia")),
    (r"cr1", ("Class_I", "LINE", "CR1")),
    (r"rte", ("Class_I", "LINE", "RTE")),
    (r"jockey|i[_-]?line|i[_-]jockey", ("Class_I", "LINE", "I-Jockey")),
    (r"\br1\b|r1[_-]", ("Class_I", "LINE", "R1")),
    (r"penelope", ("Class_I", "LINE", "Penelope")),
    (r"mariner", ("Class_II", "TIR", "TcMar-Mariner")),
    (r"\bhat\b", ("Class_II", "TIR", "hAT")),
    (r"piggybac", ("Class_II", "TIR", "PiggyBac")),
    (r"harbinger|pif", ("Class_II", "TIR", "PIF-Harbinger")),
]


def infer_te_from_dfam_model(model_name: str):
    model_lower = model_name.lower()

    for pattern, annotation in DFAM_MODEL_RULES:
        if re.search(pattern, model_lower):
            return annotation

    return ("Unknown", "Unknown", "Unknown")
