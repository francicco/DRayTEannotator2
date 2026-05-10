from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .structure_detect import TirDetection
from .tsd_from_repam import RepamTSDCall


@dataclass
class StructureEvidence:
    family_id: str
    consensus_len: int | None = None
    tir_present: bool = False
    tir_len: int = 0
    tir_identity: float = 0.0
    tir_confidence: str = "NONE"

    tsd_present: bool = False
    tsd_len: int = 0
    tsd_seq: str = "NA"
    tsd_support: float = 0.0

    structure_class: str = "NONE"
    structural_score: float = 0.0
    mite_candidate: bool = False

    mite_like_structure: bool = False
    confidence: str = "NONE"

def structure_score(
    tir_present: bool,
    tir_len: int,
    tir_identity: float,
    tsd_present: bool,
    tsd_support: float,
    consensus_len: int | None = None,
) -> float:
    score = 0.0

    if tir_present:
        score += 0.35

    if tir_len >= 15:
        score += 0.10

    if tir_len >= 25:
        score += 0.10

    if tir_identity >= 0.85:
        score += 0.15

    if tir_identity >= 0.95:
        score += 0.10

    if tsd_present:
        score += 0.25

    if tsd_support >= 0.30:
        score += 0.10

    if tsd_support >= 0.50:
        score += 0.10

    if consensus_len is not None:
        if consensus_len <= 800 and tir_present:
            score += 0.10

    return round(min(score, 1.0), 3)


def confidence_from_score(score: float) -> str:
    if score >= 0.80:
        return "HIGH"

    if score >= 0.55:
        return "MEDIUM"

    return "LOW"


def classify_structure(
    tir_present: bool,
    tsd_present: bool,
) -> str:
    if tir_present and tsd_present:
        return "TIR_TSD"

    if tir_present:
        return "TIR"

    if tsd_present:
        return "TSD"

    return "NONE"


def infer_mite_candidate(
    consensus_len: int | None,
    tir_present: bool,
    tsd_present: bool,
) -> bool:
    if consensus_len is None:
        return False

    return (
        consensus_len <= 800
        and tir_present
        and tsd_present
    )


def merge_structure_evidence(
    tir_calls: list[TirDetection],
    tsd_calls: list[RepamTSDCall],
    consensus_lengths: dict[str, int] | None = None,
) -> list[StructureEvidence]:

    if consensus_lengths is None:
        consensus_lengths = {}

    tir_map = {
        x.family_id: x
        for x in tir_calls
    }

    tsd_map = {
        x.family_id: x
        for x in tsd_calls
    }

    family_ids = sorted(
        set(tir_map) | set(tsd_map)
    )

    merged: list[StructureEvidence] = []

    for family_id in family_ids:
        tir = tir_map.get(family_id)
        tsd = tsd_map.get(family_id)

        consensus_len = consensus_lengths.get(family_id)

        tir_present = (
            tir is not None
            and tir.tir_present
        )

        tsd_present = (
            tsd is not None
            and tsd.tsd_present
        )

        score = structure_score(
            tir_present=tir_present,
            tir_len=tir.tir_len if tir else 0,
            tir_identity=tir.tir_identity if tir else 0.0,
            tsd_present=tsd_present,
            tsd_support=tsd.tsd_support if tsd else 0.0,
            consensus_len=consensus_len,
        )

        merged.append(
            StructureEvidence(
                family_id=family_id,
                consensus_len=consensus_len,
                tir_present=tir_present,
                tir_len=tir.tir_len if tir else 0,
                tir_identity=(
                    round(tir.tir_identity, 3)
                    if tir else 0.0
                ),
                tir_confidence=(
                    tir.confidence
                    if tir else "LOW"
                ),

                tsd_present=tsd_present,
                tsd_len=(
                    tsd.tsd_len
                    if tsd else 0
                ),
                tsd_seq=(
                    tsd.tsd_seq
                    if tsd and tsd.tsd_seq
                    else "NA"
                ),
                tsd_support=(
                    round(tsd.tsd_support, 3)
                    if tsd else 0.0
                ),

                structure_class=classify_structure(
                    tir_present,
                    tsd_present,
                ),

                structural_score=score,
                mite_candidate=infer_mite_candidate(
                    consensus_len=consensus_len,
                    tir_present=tir_present,
                    tsd_present=tsd_present,
                ),
                mite_like_structure=bool(tir_present and tsd_present),
                confidence=confidence_from_score(score),
            )
        )

    return merged

def write_structure_summary(
    structures: list[StructureEvidence],
    outfile: str | Path,
) -> None:

    fields = [
        "family_id",
        "consensus_len",
        "tir_present",
        "tir_len",
        "tir_identity",
        "tir_confidence",

        "tsd_present",
        "tsd_len",
        "tsd_seq",
        "tsd_support",

        "structure_class",
        "structural_score",
        "mite_candidate",
        "mite_like_structure",
        "confidence",
    ]

    with open(outfile, "w") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        for s in structures:
            row = dict(s.__dict__)

            for k, v in row.items():
                if v == "":
                    row[k] = "NA"

            writer.writerow(row)
