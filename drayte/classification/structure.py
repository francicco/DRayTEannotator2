from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from typing import Dict, Iterable


@dataclass
class StructureEvidence:
    family_id: str
    evidence_type: str
    score: float = 1.0
    source: str = "unknown"
    tsd_len: int = 0
    tsd_seq: str = ""
    tsd_support: float = 0.0
    tir_len: int = 0
    tir_identity: float = 0.0
    tir_confidence: str = "LOW"
    tir_5p_motif: str = "NA"
    tir_3p_motif: str = "NA"
    tir_terminal_motif: str = "NA"
    tir_grade: str = "ABSENT"
    structural_type: str = "unknown"
    superfamily_hint: str = "unknown"
    superfamily_confidence: str = "NONE"
    polyA_tail_len: int = 0
    polyA_fraction: float = 0.0
    poliii_a_box: bool = False
    poliii_b_box: bool = False
    poliii_score: float = 0.0
    sine_score: float = 0.0
    sine_candidate: bool = False
    sine_confidence: str = "LOW"


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


def parse_bool(v) -> bool:
    return str(v).strip().lower() in {"1", "true", "yes", "y"}


def _source_value(source: str, key: str, default: str = "") -> str:
    m = re.search(rf"(?:^|[;:]){re.escape(key)}=([^;]+)", source or "")
    if not m:
        return default
    return m.group(1)


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
                "tsd_seq": "",
                "tsd_len": 0,
                "tsd_support": 0.0,
                "tir_len": 0,
                "tir_identity": 0.0,
                "tir_confidence": "LOW",
                "tir_5p_motif": "NA",
                "tir_3p_motif": "NA",
                "tir_terminal_motif": "NA",
                "tir_grade": "ABSENT",
                "structural_type": "unknown",
                "superfamily_hint": "unknown",
                "superfamily_confidence": "NONE",
                "polyA_tail_len": 0,
                "polyA_fraction": 0.0,
                "poliii_a_box": False,
                "poliii_b_box": False,
                "poliii_score": 0.0,
                "sine_score": 0.0,
                "sine_candidate": False,
                "sine_confidence": "LOW",
            },
        )

        if ev_type == "LTR":
            fam["ltr_present"] = True
        elif ev_type == "TIR":
            fam["tir_present"] = True
            fam["tir_len"] = max(fam["tir_len"], ev.tir_len)
            fam["tir_identity"] = max(fam["tir_identity"], ev.tir_identity)
            if ev.tir_confidence in {"MEDIUM", "HIGH"}:
                fam["tir_confidence"] = ev.tir_confidence
            for motif_key in ["tir_5p_motif", "tir_3p_motif", "tir_terminal_motif"]:
                motif_value = getattr(ev, motif_key, "NA")
                if motif_value not in {"", "NA", None}:
                    fam[motif_key] = motif_value
            if ev.tir_grade not in {"", "ABSENT", "NA", None}:
                order = {"ABSENT": 0, "WEAK": 1, "MODERATE": 2, "STRONG": 3}
                if order.get(ev.tir_grade, 0) > order.get(fam["tir_grade"], 0):
                    fam["tir_grade"] = ev.tir_grade
            if ev.structural_type not in {"", "unknown", "NA", None}:
                fam["structural_type"] = ev.structural_type
            if ev.superfamily_hint not in {"", "unknown", "NA", None}:
                fam["superfamily_hint"] = ev.superfamily_hint
                fam["superfamily_confidence"] = ev.superfamily_confidence
        elif ev_type == "HELITRON":
            fam["helitron_signal"] = True
        elif ev_type == "TSD":
            fam["tsd_present"] = True
            if ev.tsd_support >= fam["tsd_support"]:
                fam["tsd_seq"] = ev.tsd_seq
                fam["tsd_len"] = ev.tsd_len
                fam["tsd_support"] = ev.tsd_support
        elif ev_type == "POLYA":
            fam["polyA_present"] = True
            fam["polyA_tail_len"] = max(fam.get("polyA_tail_len", 0), ev.polyA_tail_len)
            fam["polyA_fraction"] = max(fam.get("polyA_fraction", 0.0), ev.polyA_fraction)
            fam["poliii_a_box"] = fam.get("poliii_a_box", False) or ev.poliii_a_box
            fam["poliii_b_box"] = fam.get("poliii_b_box", False) or ev.poliii_b_box
            fam["poliii_score"] = max(fam.get("poliii_score", 0.0), ev.poliii_score)
            fam["sine_score"] = max(fam.get("sine_score", 0.0), ev.sine_score)
            fam["sine_candidate"] = fam.get("sine_candidate", False) or ev.sine_candidate
            if ev.sine_confidence in {"MEDIUM", "HIGH"}:
                fam["sine_confidence"] = ev.sine_confidence
            if ev.structural_type in {"SINE", "SINE_like"}:
                fam["structural_type"] = ev.structural_type

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
    """
    Load either the generic long structure-evidence TSV
    (family_id/evidence_type/score/source) or the structure_summary.tsv
    produced by structure_merge.py.
    """
    evidence: list[StructureEvidence] = []

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames and "evidence_type" in reader.fieldnames:
            for row in reader:
                evidence_type = normalize_structure_type(row["evidence_type"])

                if evidence_type is None:
                    continue

                source = row.get("source", "unknown")
                evidence.append(
                    StructureEvidence(
                        family_id=row["family_id"],
                        evidence_type=evidence_type,
                        score=float(row.get("score", 1.0)),
                        source=source,
                        tsd_len=int(row.get("tsd_len") or _source_value(source, "len", "0") or 0),
                        tsd_seq=row.get("tsd_seq") or _source_value(source, "seq", ""),
                        tsd_support=float(row.get("tsd_support") or 0.0),
                        tir_len=int(row.get("tir_len") or _source_value(source, "len", "0") or 0),
                        tir_identity=float(row.get("tir_identity") or _source_value(source, "identity", "0") or 0.0),
                        tir_confidence=row.get("tir_confidence") or _source_value(source, "confidence", "LOW"),
                        tir_5p_motif=row.get("tir_5p_motif") or _source_value(source, "tir_5p_motif", "NA"),
                        tir_3p_motif=row.get("tir_3p_motif") or _source_value(source, "tir_3p_motif", "NA"),
                        tir_terminal_motif=row.get("tir_terminal_motif") or _source_value(source, "tir_terminal_motif", "NA"),
                        tir_grade=row.get("tir_grade") or _source_value(source, "tir_grade", "ABSENT"),
                        structural_type=row.get("structural_type") or _source_value(source, "structural_type", "unknown"),
                        superfamily_hint=row.get("superfamily_hint") or _source_value(source, "superfamily_hint", "unknown"),
                        superfamily_confidence=row.get("superfamily_confidence") or _source_value(source, "superfamily_confidence", "NONE"),
                        polyA_tail_len=int(row.get("polyA_tail_len") or 0),
                        polyA_fraction=float(row.get("polyA_fraction") or 0.0),
                        poliii_a_box=parse_bool(row.get("poliii_a_box", False)),
                        poliii_b_box=parse_bool(row.get("poliii_b_box", False)),
                        poliii_score=float(row.get("poliii_score") or 0.0),
                        sine_score=float(row.get("sine_score") or 0.0),
                        sine_candidate=parse_bool(row.get("sine_candidate", False)),
                        sine_confidence=row.get("sine_confidence") or "LOW",
                    )
                )
            return evidence

        for row in reader:
            family_id = row["family_id"]

            if parse_bool(row.get("tir_present", False)):
                evidence.append(
                    StructureEvidence(
                        family_id=family_id,
                        evidence_type="TIR",
                        score=float(row.get("tir_identity") or 1.0),
                        source="structure_summary:TIR",
                        tir_len=int(row.get("tir_len") or 0),
                        tir_identity=float(row.get("tir_identity") or 0.0),
                        tir_confidence=row.get("tir_confidence") or "LOW",
                        tir_5p_motif=row.get("tir_5p_motif") or "NA",
                        tir_3p_motif=row.get("tir_3p_motif") or "NA",
                        tir_terminal_motif=row.get("tir_terminal_motif") or "NA",
                        tir_grade=row.get("tir_grade") or "ABSENT",
                        structural_type=row.get("structural_type") or "unknown",
                        superfamily_hint=row.get("superfamily_hint") or "unknown",
                        superfamily_confidence=row.get("superfamily_confidence") or "NONE",
                    )
                )

            if parse_bool(row.get("tsd_present", False)):
                evidence.append(
                    StructureEvidence(
                        family_id=family_id,
                        evidence_type="TSD",
                        score=float(row.get("tsd_support") or 1.0),
                        source="structure_summary:TSD",
                        tsd_len=int(row.get("tsd_len") or 0),
                        tsd_seq=(row.get("tsd_seq") or "").replace("NA", ""),
                        tsd_support=float(row.get("tsd_support") or 0.0),
                    )
                )

            if parse_bool(row.get("polyA_present", False)) or parse_bool(row.get("sine_candidate", False)):
                evidence.append(
                    StructureEvidence(
                        family_id=family_id,
                        evidence_type="POLYA",
                        score=float(row.get("sine_score") or row.get("polyA_fraction") or 1.0),
                        source="structure_summary:SINE",
                        structural_type=row.get("structural_type") or "unknown",
                        polyA_tail_len=int(row.get("polyA_tail_len") or 0),
                        polyA_fraction=float(row.get("polyA_fraction") or 0.0),
                        poliii_a_box=parse_bool(row.get("poliii_a_box", False)),
                        poliii_b_box=parse_bool(row.get("poliii_b_box", False)),
                        poliii_score=float(row.get("poliii_score") or 0.0),
                        sine_score=float(row.get("sine_score") or 0.0),
                        sine_candidate=parse_bool(row.get("sine_candidate", False)),
                        sine_confidence=row.get("sine_confidence") or "LOW",
                    )
                )

    return evidence


def write_structure_evidence_tsv(evidence, outpath):
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

