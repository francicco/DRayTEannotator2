from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ParsedClassification:
    family_id: str
    rm_class: str
    rm_superfamily: str


ORDER_MAP = {
    "LTR": "LTR",
    "LINE": "LINE",
    "SINE": "SINE",
    "DNA": "TIR",
    "RC": "Helitron",
    "Helitron": "Helitron",
}


CLASS_MAP = {
    "LTR": "Class_I",
    "LINE": "Class_I",
    "SINE": "Class_I",
    "TIR": "Class_II",
    "Helitron": "Class_II",
}


def parse_repeatmasker_header(header: str) -> ParsedClassification | None:

    if "#" not in header:
        return None

    family_id, annotation = header.split("#", 1)

    if "/" in annotation:
        rm_class, rm_superfamily = annotation.split("/", 1)
    else:
        rm_class = annotation
        rm_superfamily = "Unknown"

    return ParsedClassification(
        family_id=family_id,
        rm_class=rm_class,
        rm_superfamily=rm_superfamily,
    )


def classification_from_repeatmasker(parsed: ParsedClassification):

    order = ORDER_MAP.get(
        parsed.rm_class,
        parsed.rm_class,
    )

    te_class = CLASS_MAP.get(
        order,
        "Unknown",
    )

    return {
        "family_id": parsed.family_id,
        "expected_class": te_class,
        "expected_order": order,
        "expected_superfamily": parsed.rm_superfamily,
    }
