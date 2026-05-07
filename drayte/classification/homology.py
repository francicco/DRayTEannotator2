from __future__ import annotations

from dataclasses import dataclass

from .parsers.repeatmasker import parse_repeatmasker_header


@dataclass
class HomologyEvidence:
    family_id: str
    homology_class: str
    homology_order: str
    homology_superfamily: str
    homology_score: float


def infer_te_class_from_rm_class(rm_class: str) -> str:
    if rm_class in {"LINE", "LTR", "SINE"}:
        return "Class_I"

    if rm_class in {"DNA", "RC", "Helitron"}:
        return "Class_II"

    return "Unknown"


def infer_order_from_rm_class(rm_class: str) -> str:
    if rm_class == "DNA":
        return "TIR"

    if rm_class in {"RC", "Helitron"}:
        return "Helitron"

    if rm_class in {"LINE", "LTR", "SINE"}:
        return rm_class

    return "Unknown"


def homology_from_repeatmasker_header(
    seq_id: str,
    score: float = 1.0,
) -> HomologyEvidence | None:

    parsed = parse_repeatmasker_header(seq_id)

    if parsed is None:
        return None

    return HomologyEvidence(
        family_id=parsed.family_id,
        homology_class=infer_te_class_from_rm_class(parsed.rm_class),
        homology_order=infer_order_from_rm_class(parsed.rm_class),
        homology_superfamily=parsed.rm_superfamily,
        homology_score=score,
    )
