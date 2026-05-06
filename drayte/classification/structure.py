from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable

@dataclass

class StructureEvidence:
    family_id: str
    evidence_type: str
    score: float = 1.0
    source: str = "unknown"

VALID_STRUCTURE_TYPES = {
    "LTR",
    "TIR",
    "HELITRON",
    "TSD",
    "POLYA",
}

def normalize_structure_type(value: str) -> str | None:
    v = value.strip().upper()
    if v in {"LTR", "LTR_RETROTRANSPOSON", "LTR_STRUCTURE"}:
        return "LTR"
    if v in {"TIR", "TIR_STRUCTURE", "DNA_TIR"}:
        return "TIR"
    if v in {"HELITRON", "RC", "ROLLING_CIRCLE", "HELITRON_SIGNAL"}:
        return "HELITRON"
    if v in {"TSD", "TARGET_SITE_DUPLICATION"}:
        return "TSD"
    if v in {"POLYA", "POLY_A", "POLYA_TAIL"}:
        return "POLYA"
    return None

def summarize_structure_evidence(
    evidence: Iterable[StructureEvidence],
) -> Dict[str, dict]:
    summary: Dict[str, dict] = {}
    for ev in evidence:
        ev_type = normalize_structure_type(ev.evidence_type)
        if ev_type is None:
            continue
        fam = summary.setdefault(
            ev.family_id,
            {
                "ltr_present": False,
                "tir_present": False,
                "helitron_signal": False,
                "tsd_present": False,
                "polyA_present": False,
                "structure_sources": set(),
                "best_structure_score": 0.0,
            },
        )

        if ev_type == "LTR":
            fam["ltr_present"] = True
        elif ev_type == "TIR":
            fam["tir_present"] = True
        elif ev_type == "HELITRON":
            fam["helitron_signal"] = True
        elif ev_type == "TSD":
            fam["tsd_present"] = True
        elif ev_type == "POLYA":
            fam["polyA_present"] = True
        fam["structure_sources"].add(ev.source)
        fam["best_structure_score"] = max(fam["best_structure_score"], ev.score)
    for fam in summary.values():
        fam["structure_sources"] = sorted(fam["structure_sources"])
    return summary


def evidence_from_structure_candidates(candidates) -> list[StructureEvidence]:
    evidence: list[StructureEvidence] = []

    for cand in candidates:
        ev_type = normalize_structure_type(cand.module)

        if ev_type is None:
            continue

        evidence.append(
            StructureEvidence(
                family_id=cand.candidate_id,
                evidence_type=ev_type,
                score=1.0,
                source=f"structure:{cand.module}",
            )
        )

    return evidence


def load_structure_evidence_tsv(path) -> list[StructureEvidence]:
    import csv

    evidence = []

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            evidence_type = normalize_structure_type(
                row["evidence_type"]
            )

            if evidence_type is None:
                continue

            evidence.append(
                StructureEvidence(
                    family_id=row["family_id"],
                    evidence_type=evidence_type,
                    score=float(row.get("score", 1.0)),
                    source=row.get("source", "unknown"),
                )
            )

    return evidence


def write_structure_evidence_tsv(evidence, outpath):
    import csv

    fields = [
        "family_id",
        "evidence_type",
        "score",
        "source",
    ]

    with open(outpath, "w") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        for ev in evidence:
            writer.writerow({
                "family_id": ev.family_id,
                "evidence_type": ev.evidence_type,
                "score": ev.score,
                "source": ev.source,
            })
