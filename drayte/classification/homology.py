from __future__ import annotations

from dataclasses import dataclass

from .parsers.repeatmasker import parse_repeatmasker_header


@dataclass
class HomologyEvidence:
    family_id: str
    homology_class: str
    homology_superfamily: str
    homology_score: float


CLASS_NORMALIZATION = {
    "RC": "Helitron",
}


def normalize_homology_class(name: str) -> str:
    return CLASS_NORMALIZATION.get(name, name)


def homology_from_repeatmasker_header(
    seq_id: str,
    score: float = 1.0,
) -> HomologyEvidence | None:

    parsed = parse_repeatmasker_header(seq_id)

    if parsed is None:
        return None

    return HomologyEvidence(
        family_id=parsed.family_id,
        homology_class=normalize_homology_class(
            parsed.rm_class
        ),
        homology_superfamily=parsed.rm_superfamily,
        homology_score=score,
    )
