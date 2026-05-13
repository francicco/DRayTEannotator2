from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from drayte.classification.complexity_detect import ComplexityDetection
from drayte.classification.sine_detect import SineDetection
from drayte.classification.structure_detect import TirDetection
from drayte.classification.tsd_from_repam import RepamTSDCall
from .ltr_detect import LtrDetection

@dataclass
class StructureEvidence:
    family_id: str
    consensus_len: int | None = None

    seq_entropy: float = 0.0
    dominant_1mer_fraction: float = 0.0
    dominant_2mer_fraction: float = 0.0
    dominant_3mer_fraction: float = 0.0
    dominant_4mer_fraction: float = 0.0
    max_homopolymer_run: int = 0
    low_complexity_candidate: bool = False
    simple_repeat_candidate: bool = False
    min_window_entropy: float = 0.0
    max_window_dom1_fraction: float = 0.0
    max_window_homopolymer_run: int = 0

    is_low_complexity: bool = False
    is_simple_repeat: bool = False
    complexity_score: float = 0.0
    complexity_class: str = "none"
    complexity_confidence: str = "LOW"

    tir_present: bool = False
    tir_len: int = 0
    tir_identity: float = 0.0
    tir_confidence: str = "NONE"
    tir_5p_motif: str = "NA"
    tir_3p_motif: str = "NA"
    tir_terminal_motif: str = "NA"

    tsd_present: bool = False
    tsd_len: int = 0
    tsd_seq: str = "NA"
    tsd_support: float = 0.0
    tsd_len_support: float = 0.0
    tsd_exact_support: float = 0.0
    tsd_confidence: str = "NONE"

    polyA_present: bool = False
    polyA_tail_len: int = 0
    polyA_fraction: float = 0.0
    tail_type: str = "none"

    poliii_a_box: bool = False
    poliii_b_box: bool = False
    poliii_score: float = 0.0

    sine_head_type: str = "unknown"
    sine_head_score: float = 0.0
    sine_score: float = 0.0
    sine_candidate: bool = False
    sine_confidence: str = "LOW"

    ltr_present: bool = False
    ltr_len: int = 0
    ltr_identity: float = 0.0
    ltr_internal_len: int = 0
    ltr_score: float = 0.0
    ltr_structural_type: str = "none"
    ltr_confidence: str = "LOW"

    tg_ca_motif: bool = False
    left_ltr_motif: str = "NA"
    right_ltr_motif: str = "NA"
    pbs_like: bool = False
    ppt_like: bool = False

    structure_class: str = "NONE"
    structural_type: str = "unknown"
    tir_grade: str = "ABSENT"
    superfamily_hint: str = "unknown"
    superfamily_confidence: str = "NONE"
    structural_score: float = 0.0
    mite_candidate: bool = False
    tir_tsd_element: bool = False
    mite_like_structure: bool = False
    confidence: str = "NONE"


def structure_score(tir_present, tir_len, tir_identity, tsd_present, tsd_support, consensus_len=None) -> float:
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
    if consensus_len is not None and consensus_len <= 800 and tir_present:
        score += 0.10
    return round(min(score, 1.0), 3)


def confidence_from_score(score: float) -> str:
    if score >= 0.80:
        return "HIGH"
    if score >= 0.55:
        return "MEDIUM"
    return "LOW"


def classify_structure(tir_present: bool, tsd_present: bool) -> str:
    if tir_present and tsd_present:
        return "TIR_TSD"
    if tir_present:
        return "TIR"
    if tsd_present:
        return "TSD"
    return "NONE"


def infer_tir_grade(tir_present, tir_len, tir_identity, tir_confidence="LOW") -> str:
    if not tir_present:
        return "ABSENT"
    conf = (tir_confidence or "LOW").upper()
    if conf == "HIGH" and tir_len >= 20 and tir_identity >= 0.90:
        return "STRONG"
    if conf in {"HIGH", "MEDIUM"} and tir_len >= 15 and tir_identity >= 0.85:
        return "MODERATE"
    return "WEAK"


def infer_mite_candidate(consensus_len, tir_present, tsd_present, tir_grade="ABSENT") -> bool:
    return bool(
        consensus_len is not None
        and consensus_len <= 800
        and tir_present
        and tsd_present
        and tir_grade in {"MODERATE", "STRONG"}
    )


def infer_structural_type(
    consensus_len: int | None,
    tir_present: bool,
    tsd_present: bool,
    tir_grade: str,
    sine_candidate: bool = False,
    sine_score: float = 0.0,
    ltr_present: bool = False,
    ltr_structural_type: str = "none",
    simple_repeat_candidate: bool = False,
    low_complexity_candidate: bool = False,
) -> str:
    is_short = consensus_len is not None and consensus_len <= 800
    has_good_tir = tir_grade in {"MODERATE", "STRONG"}

    if tir_present:
        if tsd_present:
            if is_short and has_good_tir:
                return "MITE"
            if is_short:
                return "short_TIR_TSD_element"
            return "TIR_TSD_element"
        if is_short:
            return "short_TIR_element"
        return "TIR_element"

    if sine_candidate and not tsd_present:
        return "SINE"

    if (
        not tsd_present
        and consensus_len is not None
        and consensus_len <= 700
        and sine_score >= 0.50
    ):
        return "SINE_like"

    if ltr_present and ltr_structural_type != "none":
        return ltr_structural_type

    if simple_repeat_candidate:
        return "simple_repeat"

    if low_complexity_candidate:
        return "low_complexity"

    return "unknown"


def infer_tir_tsd_element(tir_present: bool, tsd_present: bool) -> bool:
    return bool(tir_present and tsd_present)


def _clean_motif(seq: str, motif_len: int = 8) -> str:
    seq = (seq or "").upper()
    seq = "".join(base for base in seq if base in {"A", "C", "G", "T"})
    return seq[:motif_len] if seq else "NA"


def _revcomp(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return (seq or "").translate(table)[::-1].upper()


def extract_tir_terminal_motifs(tir: TirDetection | None, motif_len: int = 8) -> tuple[str, str, str]:
    if tir is None or not tir.tir_present:
        return "NA", "NA", "NA"

    left = _clean_motif(tir.left_tir_seq, motif_len)
    right = _clean_motif(_revcomp(tir.right_tir_seq), motif_len)

    if left != "NA" and right != "NA":
        shared = []
        for a, b in zip(left, right):
            if a == b:
                shared.append(a)
            else:
                break
        terminal = "".join(shared) if shared else left
    else:
        terminal = left if left != "NA" else right

    return left, right, terminal


def motif_suggests_cacta(motif: str) -> bool:
    motif = (motif or "").upper()
    return motif.startswith(("CACTA", "CACTG", "CCC"))


def motif_suggests_tc1_mariner(motif: str) -> bool:
    motif = (motif or "").upper()
    return motif.startswith(("TA", "CA", "TG"))


def motif_suggests_piggybac(motif: str) -> bool:
    motif = (motif or "").upper()
    return motif.startswith(("TTAA", "TTA", "CCA"))


def infer_tir_superfamily_hint(
    tsd_present: bool,
    tsd_len: int,
    tsd_seq: str,
    tir_present: bool = True,
    tir_grade: str = "ABSENT",
    tir_terminal_motif: str = "NA",
) -> tuple[str, str]:
    if not tir_present or not tsd_present or tsd_len <= 0:
        return "unknown", "NONE"

    seq = (tsd_seq or "").upper()
    motif = (tir_terminal_motif or "NA").upper()
    grade = (tir_grade or "ABSENT").upper()
    strong_tir = grade in {"MODERATE", "STRONG"}

    if tsd_len == 2 and seq == "TA":
        return "Tc1-Mariner_like", "HIGH" if strong_tir or motif_suggests_tc1_mariner(motif) else "MEDIUM"

    if tsd_len == 4 and seq == "TTAA":
        return "piggyBac_like", "HIGH" if motif_suggests_piggybac(motif) or strong_tir else "MEDIUM"

    if tsd_len == 3:
        if seq in {"TAA", "TTA"}:
            return "PIF-Harbinger_like", "MEDIUM" if strong_tir else "LOW"
        if motif_suggests_cacta(motif):
            return "CACTA_like", "MEDIUM"
        return "PIF-Harbinger_or_CACTA_like", "LOW"

    if tsd_len == 2:
        if motif_suggests_cacta(motif):
            return "CACTA_like", "MEDIUM" if strong_tir else "LOW"
        return "2bp_TSD_TIR_ambiguous", "LOW"

    if tsd_len == 5:
        return "Transib_like", "LOW"
    if tsd_len == 6:
        return "Maverick_or_LTR_like", "LOW"
    if tsd_len == 8:
        return "hAT_P_Merlin_ambiguous", "LOW"
    if tsd_len in {9, 10, 11}:
        return "Mutator_like", "MEDIUM" if strong_tir else "LOW"

    return "unknown", "NONE"


def merge_structure_evidence(
    tir_calls: list[TirDetection],
    tsd_calls: list[RepamTSDCall],
    consensus_lengths: dict[str, int] | None = None,
    complexity_calls: list[ComplexityDetection] | None = None,
    sine_calls: list[SineDetection] | None = None,
    ltr_calls: list[LtrDetection] | None = None,
) -> list[StructureEvidence]:
    if consensus_lengths is None:
        consensus_lengths = {}

    tir_map = {x.family_id: x for x in tir_calls}
    tsd_map = {x.family_id: x for x in tsd_calls}
    sine_map = {x.family_id: x for x in (sine_calls or [])}
    complexity_map = {x.family_id: x for x in (complexity_calls or [])}
    ltr_map = {x.family_id: x for x in (ltr_calls or [])}

    family_ids = sorted(set(tir_map) | set(tsd_map) | set(sine_map) | set(complexity_map) | set(ltr_map))
    merged: list[StructureEvidence] = []

    for family_id in family_ids:
        tir = tir_map.get(family_id)
        tsd = tsd_map.get(family_id)
        sine = sine_map.get(family_id)
        ltr = ltr_map.get(family_id)
        complexity = complexity_map.get(family_id)

        consensus_len = consensus_lengths.get(family_id)
        if consensus_len is None and sine is not None and sine.consensus_len:
            consensus_len = sine.consensus_len
        if consensus_len is None and complexity is not None and complexity.consensus_len:
            consensus_len = complexity.consensus_len

        tir_present = bool(tir and tir.tir_present)
        tsd_present = bool(tsd and tsd.tsd_present)

        tir_len = tir.tir_len if tir else 0
        tir_identity = tir.tir_identity if tir else 0.0
        tir_confidence = tir.confidence if tir else "LOW"
        tir_5p_motif, tir_3p_motif, tir_terminal_motif = extract_tir_terminal_motifs(tir)

        tsd_len = tsd.tsd_len if tsd else 0
        tsd_seq = tsd.tsd_seq if tsd and tsd.tsd_seq else "NA"
        tsd_support = tsd.tsd_support if tsd else 0.0
        tsd_len_support = getattr(tsd, "tsd_len_support", tsd_support) if tsd else 0.0
        tsd_exact_support = getattr(tsd, "tsd_exact_support", tsd_support) if tsd else 0.0
        tsd_confidence = getattr(tsd, "tsd_confidence", "NONE") if tsd else "NONE"

        polyA_present = sine.polyA_present if sine else False
        polyA_tail_len = sine.polyA_tail_len if sine else 0
        polyA_fraction = sine.polyA_fraction if sine else 0.0
        tail_type = sine.tail_type if sine else "none"
        poliii_a_box = sine.poliii_a_box if sine else False
        poliii_b_box = sine.poliii_b_box if sine else False
        poliii_score = sine.poliii_score if sine else 0.0
        sine_head_type = sine.sine_head_type if sine else "unknown"
        sine_head_score = sine.sine_head_score if sine else 0.0
        raw_sine_score = sine.sine_score if sine else 0.0
        raw_sine_candidate = sine.sine_candidate if sine else False
        sine_confidence = sine.confidence if sine else "LOW"

        ltr_present = ltr.ltr_present if ltr else False
        ltr_len = ltr.ltr_len if ltr else 0
        ltr_identity = ltr.ltr_identity if ltr else 0.0
        ltr_internal_len = ltr.internal_len if ltr else 0
        ltr_score = ltr.ltr_score if ltr else 0.0
        ltr_structural_type = ltr.ltr_structural_type if ltr else "none"
        ltr_confidence = ltr.confidence if ltr else "LOW"

        tg_ca_motif = ltr.tg_ca_motif if ltr else False
        left_ltr_motif = ltr.left_ltr_motif if ltr else "NA"
        right_ltr_motif = ltr.right_ltr_motif if ltr else "NA"
        pbs_like = ltr.pbs_like if ltr else False
        ppt_like = ltr.ppt_like if ltr else False

        seq_entropy = complexity.seq_entropy if complexity else 0.0
        dominant_1mer_fraction = complexity.dominant_1mer_fraction if complexity else 0.0
        dominant_2mer_fraction = complexity.dominant_2mer_fraction if complexity else 0.0
        dominant_3mer_fraction = complexity.dominant_3mer_fraction if complexity else 0.0
        dominant_4mer_fraction = complexity.dominant_4mer_fraction if complexity else 0.0
        max_homopolymer_run = complexity.max_homopolymer_run if complexity else 0
        low_complexity_candidate = complexity.low_complexity_candidate if complexity else False
        simple_repeat_candidate = complexity.simple_repeat_candidate if complexity else False
        complexity_class = complexity.complexity_class if complexity else "none"
        complexity_confidence = complexity.complexity_confidence if complexity else "LOW"

        min_win_entropy = (
            complexity.min_window_entropy
            if complexity and hasattr(complexity, "min_window_entropy")
            else 0.0
        )
        
        max_win_dom1 = (
            complexity.max_window_dom1_fraction
            if complexity and hasattr(complexity, "max_window_dom1_fraction")
            else 0.0
        )
        
        max_win_homo = (
            complexity.max_window_homopolymer_run
            if complexity and hasattr(complexity, "max_window_homopolymer_run")
            else 0
        )

        tir_grade = infer_tir_grade(tir_present, tir_len, tir_identity, tir_confidence)

        sine_candidate = bool(raw_sine_candidate and not tir_present and not tsd_present)
        sine_score = 0.0 if tir_present else raw_sine_score

        structural_type = infer_structural_type(
            consensus_len=consensus_len,
            tir_present=tir_present,
            tsd_present=tsd_present,
            tir_grade=tir_grade,
            sine_candidate=sine_candidate,
            sine_score=sine_score,
            simple_repeat_candidate=simple_repeat_candidate,
            low_complexity_candidate=low_complexity_candidate,
            ltr_present=ltr_present,
            ltr_structural_type=ltr_structural_type,
        )

        superfamily_hint, superfamily_confidence = infer_tir_superfamily_hint(
            tsd_present=tsd_present,
            tsd_len=tsd_len,
            tsd_seq=tsd_seq,
            tir_present=tir_present,
            tir_grade=tir_grade,
            tir_terminal_motif=tir_terminal_motif,
        )

        score = structure_score(
            tir_present=tir_present,
            tir_len=tir_len,
            tir_identity=tir_identity,
            tsd_present=tsd_present,
            tsd_support=tsd_support,
            consensus_len=consensus_len,
        )

        merged.append(
            StructureEvidence(
                family_id=family_id,
                consensus_len=consensus_len,
                seq_entropy=seq_entropy,
                dominant_1mer_fraction=dominant_1mer_fraction,
                dominant_2mer_fraction=dominant_2mer_fraction,
                dominant_3mer_fraction=dominant_3mer_fraction,
                dominant_4mer_fraction=dominant_4mer_fraction,
                max_homopolymer_run=max_homopolymer_run,
                low_complexity_candidate=low_complexity_candidate,
                min_window_entropy=min_win_entropy,
                max_window_dom1_fraction=max_win_dom1,
                max_window_homopolymer_run=max_win_homo,
                simple_repeat_candidate=simple_repeat_candidate,
                complexity_class=complexity_class,
                complexity_confidence=complexity_confidence,
                tir_present=tir_present,
                tir_len=tir_len,
                tir_identity=round(tir_identity, 3),
                tir_confidence=tir_confidence,
                tir_5p_motif=tir_5p_motif,
                tir_3p_motif=tir_3p_motif,
                tir_terminal_motif=tir_terminal_motif,
                tsd_present=tsd_present,
                tsd_len=tsd_len,
                tsd_seq=tsd_seq,
                tsd_support=round(tsd_support, 3),
                tsd_len_support=round(tsd_len_support, 3),
                tsd_exact_support=round(tsd_exact_support, 3),
                tsd_confidence=tsd_confidence,
                polyA_present=polyA_present,
                polyA_tail_len=polyA_tail_len,
                polyA_fraction=round(polyA_fraction, 3),
                tail_type=tail_type,
                poliii_a_box=poliii_a_box,
                poliii_b_box=poliii_b_box,
                poliii_score=round(poliii_score, 3),
                sine_head_type=sine_head_type,
                sine_head_score=round(sine_head_score, 3),
                sine_score=round(sine_score, 3),
                sine_candidate=sine_candidate,
                sine_confidence=sine_confidence,
                ltr_present=ltr_present,
                ltr_len=ltr_len,
                ltr_identity=round(ltr_identity, 3),
                ltr_internal_len=ltr_internal_len,
                ltr_score=round(ltr_score, 3),
                ltr_structural_type=ltr_structural_type,
                ltr_confidence=ltr_confidence,
                tg_ca_motif=tg_ca_motif,
                left_ltr_motif=left_ltr_motif,
                right_ltr_motif=right_ltr_motif,
                pbs_like=pbs_like,
                ppt_like=ppt_like,
                structure_class=classify_structure(tir_present, tsd_present),
                structural_type=structural_type,
                tir_grade=tir_grade,
                superfamily_hint=superfamily_hint,
                superfamily_confidence=superfamily_confidence,
                structural_score=score,
                mite_candidate=infer_mite_candidate(consensus_len, tir_present, tsd_present, tir_grade),
                tir_tsd_element=infer_tir_tsd_element(tir_present, tsd_present),
                mite_like_structure=infer_tir_tsd_element(tir_present, tsd_present),
                confidence=confidence_from_score(score),
            )
        )

    return merged


def write_structure_summary(structures: list[StructureEvidence], outfile: str | Path) -> None:
    fields = [
        "family_id",
        "consensus_len",
        "seq_entropy",
        "dominant_1mer_fraction",
        "dominant_2mer_fraction",
        "dominant_3mer_fraction",
        "dominant_4mer_fraction",
        "max_homopolymer_run",
        "low_complexity_candidate",
        "simple_repeat_candidate",
        "complexity_class",
        "min_window_entropy",
        "max_window_dom1_fraction",
        "max_window_homopolymer_run",
        "complexity_confidence",
        "tir_present",
        "tir_len",
        "tir_identity",
        "tir_confidence",
        "tir_5p_motif",
        "tir_3p_motif",
        "tir_terminal_motif",
        "tsd_present",
        "tsd_len",
        "tsd_seq",
        "tsd_support",
        "tsd_len_support",
        "tsd_exact_support",
        "tsd_confidence",
        "polyA_present",
        "polyA_tail_len",
        "polyA_fraction",
        "tail_type",
        "poliii_a_box",
        "poliii_b_box",
        "poliii_score",
        "sine_head_type",
        "sine_head_score",
        "sine_score",
        "sine_candidate",
        "sine_confidence",
        "ltr_present",
        "ltr_len",
        "ltr_identity",
        "ltr_internal_len",
        "ltr_score",
        "ltr_structural_type",
        "ltr_confidence",
        "tg_ca_motif",
        "left_ltr_motif",
        "right_ltr_motif",
        "pbs_like",
        "ppt_like",
        "structure_class",
        "structural_type",
        "tir_grade",
        "superfamily_hint",
        "superfamily_confidence",
        "structural_score",
        "mite_candidate",
        "tir_tsd_element",
        "mite_like_structure",
        "confidence",
    ]

    with open(outfile, "w") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for s in structures:
            row = dict(s.__dict__)
            for k, v in row.items():
                if v == "":
                    row[k] = "NA"
            writer.writerow(row)
