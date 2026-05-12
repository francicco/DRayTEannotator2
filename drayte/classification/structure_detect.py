from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from .ids import clean_family_id
from .structure_refine import refine_tir_detection

@dataclass
class TirDetection:
    family_id: str

    tir_present: bool = False
    tir_len: int = 0
    tir_identity: float = 0.0

    left_tir_start: int = 0
    left_tir_end: int = 0
    right_tir_start: int = 0
    right_tir_end: int = 0

    left_tir_seq: str = ""
    right_tir_seq: str = ""

    tsd_present: bool = False
    tsd_len: int = 0
    tsd_seq: str = ""

    terminality_score: float = 0.0
    complexity_score: float = 0.0

    window: int = 0
    confidence: str = "LOW"


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def is_simple_repeat(seq: str) -> bool:
    seq = seq.upper()

    for k in [1, 2, 3]:
        if len(seq) < k * 3:
            continue

        motif = seq[:k]

        rebuilt = (motif * ((len(seq) // k) + 1))[:len(seq)]

        if rebuilt == seq:
            return True

    return False

def sequence_complexity(seq: str) -> float:
    seq = seq.upper()
    seq = "".join(x for x in seq if x in {"A", "C", "G", "T"})

    if not seq:
        return 0.0

    return len(set(seq)) / 4.0


def is_low_complexity(
    seq: str,
    min_complexity: float = 0.50,
    max_single_base_fraction: float = 0.75,
) -> bool:
    seq = seq.upper()
    seq = "".join(x for x in seq if x in {"A", "C", "G", "T"})

    if not seq:
        return True

    complexity = sequence_complexity(seq)

    max_base_fraction = max(
        seq.count(base) / len(seq)
        for base in "ACGT"
    )

    return (
        complexity < min_complexity
        or max_base_fraction > max_single_base_fraction
    )


def pair_identity(a: str, b: str) -> float:
    if len(a) != len(b) or not a:
        return 0.0

    matches = sum(
        1 for x, y in zip(a.upper(), b.upper())
        if x == y
    )

    return matches / len(a)


def detect_tsd_near_tir(
    seq: str,
    left_tir_start: int,
    right_tir_end: int,
    min_len: int = 2,
    max_len: int = 10,
) -> tuple[bool, int, str]:
    """
    Detect TSD outside the inferred TIRs.

    Expected structure:

        left_flank | TSD | left_TIR ... right_TIR | TSD | right_flank

    Coordinates are zero-based internally.
    """
    seq = seq.upper()

    for k in range(max_len, min_len - 1, -1):
        if left_tir_start < k:
            continue

        if right_tir_end + k > len(seq):
            continue

        left_tsd = seq[left_tir_start - k:left_tir_start]
        right_tsd = seq[right_tir_end:right_tir_end + k]

        if left_tsd == right_tsd and not is_low_complexity(left_tsd):
            return True, k, left_tsd

    return False, 0, ""


def candidate_is_better(candidate: dict, best: dict | None) -> bool:
    """
    Ranking strategy for TIR candidates.

    We prefer biologically stronger evidence over simply longer matches:
    TSD support > identity > terminality > TIR length > smaller window.
    """
    if best is None:
        return True

    if candidate.get("tsd_present", False) != best.get("tsd_present", False):
        return candidate.get("tsd_present", False)

    if candidate["tir_identity"] != best["tir_identity"]:
        return candidate["tir_identity"] > best["tir_identity"]

    if candidate["terminality_score"] != best["terminality_score"]:
        return candidate["terminality_score"] > best["terminality_score"]

    if candidate["tir_len"] != best["tir_len"]:
        return candidate["tir_len"] > best["tir_len"]

    return candidate.get("window", 10**9) < best.get("window", 10**9)


def find_best_terminal_tir(
    seq: str,
    window: int = 120,
    min_len: int = 16,
    max_len: int = 60,
    min_identity: float = 0.90,
    min_complexity: float = 0.50,
) -> dict | None:
    """
    Search only terminal windows.

    The left TIR is searched in the 5' terminal window.
    The right TIR is searched in the reverse complement of the 3' terminal
    window. This avoids whole-sequence self-similarity overcalling.
    """
    seq = seq.upper()
    n = len(seq)

    if n < min_len * 2:
        return None

    terminal_window = min(window, n // 2)

    left_window = seq[:terminal_window]
    right_window_raw = seq[n - terminal_window:]
    right_window_rc = revcomp(right_window_raw)

    best = None

    for i in range(0, len(left_window) - min_len + 1):
        for j in range(0, len(right_window_rc) - min_len + 1):
            max_k = min(
                len(left_window) - i,
                len(right_window_rc) - j,
                max_len,
            )

            for k in range(min_len, max_k + 1):
                left_tir = left_window[i:i + k]
                right_tir_rc = right_window_rc[j:j + k]

                if is_simple_repeat(left_tir):
                    continue
                
                if is_simple_repeat(right_tir_rc):
                    continue

                if is_low_complexity(
                    left_tir,
                    min_complexity=min_complexity,
                ):
                    continue

                if is_low_complexity(
                    right_tir_rc,
                    min_complexity=min_complexity,
                ):
                    continue

                identity = pair_identity(left_tir, right_tir_rc)

                if identity < min_identity:
                    continue

                # Convert coordinates from right-window reverse-complement space
                # back to original sequence coordinates.
                right_tir_start = n - terminal_window + (
                    terminal_window - j - k
                )
                right_tir_end = right_tir_start + k

                left_tir_start = i
                left_tir_end = i + k

                left_terminality = 1.0 - (
                    left_tir_start / terminal_window
                )
                right_terminality = right_tir_end / n
                terminality_score = (
                    left_terminality + right_terminality
                ) / 2.0

                if left_tir_start > 10:
                    continue
                
                if (n - right_tir_end) > 10:
                    continue

                complexity_score = min(
                    sequence_complexity(left_tir),
                    sequence_complexity(right_tir_rc),
                )

                tsd_present, tsd_len, tsd_seq = detect_tsd_near_tir(
                    seq=seq,
                    left_tir_start=left_tir_start,
                    right_tir_end=right_tir_end,
                )

                candidate = {
                    "tir_len": k,
                    "tir_identity": identity,
                    "left_tir_start": left_tir_start,
                    "left_tir_end": left_tir_end,
                    "right_tir_start": right_tir_start,
                    "right_tir_end": right_tir_end,
                    "left_tir_seq": left_tir,
                    "right_tir_seq": seq[right_tir_start:right_tir_end],
                    "terminality_score": terminality_score,
                    "complexity_score": complexity_score,
                    "tsd_present": tsd_present,
                    "tsd_len": tsd_len,
                    "tsd_seq": tsd_seq,
                    "window": window,
                }

                if candidate_is_better(candidate, best):
                    best = candidate

    return best


def confidence_from_tir(
    tir_len: int,
    tir_identity: float,
    tsd_present: bool,
    terminality_score: float,
    complexity_score: float,
) -> str:
    score = 0.0

    if tir_len >= 18:
        score += 0.20

    if tir_len >= 30:
        score += 0.10

    if tir_identity >= 0.85:
        score += 0.25

    if tir_identity >= 0.95:
        score += 0.10

    if tsd_present:
        score += 0.25

    if terminality_score >= 0.90:
        score += 0.10

    if complexity_score >= 0.50:
        score += 0.10

    if score >= 0.80:
        return "HIGH"

    if score >= 0.55:
        return "MEDIUM"

    return "LOW"


def best_terminal_inverted_repeat(
    seq: str,
    family_id: str = "unknown",
    window: int = 120,
    min_len: int = 12,
    max_len: int = 60,
    min_identity: float = 0.85,
    min_complexity: float = 0.50,
    logger=None,
) -> TirDetection | None:
    seq = seq.upper()

    best = find_best_terminal_tir(
        seq=seq,
        window=window,
        min_len=min_len,
        max_len=max_len,
        min_identity=min_identity,
        min_complexity=min_complexity,
    )

    if best is None:
        return None

    confidence = confidence_from_tir(
        tir_len=best["tir_len"],
        tir_identity=best["tir_identity"],
        tsd_present=best["tsd_present"],
        terminality_score=best["terminality_score"],
        complexity_score=best["complexity_score"],
    )

    det = TirDetection(
        family_id=family_id,
        tir_present=True,
        tir_len=best["tir_len"],
        tir_identity=best["tir_identity"],
        left_tir_start=best["left_tir_start"] + 1,
        left_tir_end=best["left_tir_end"],
        right_tir_start=best["right_tir_start"] + 1,
        right_tir_end=best["right_tir_end"],
        left_tir_seq=best["left_tir_seq"],
        right_tir_seq=best["right_tir_seq"],
        tsd_present=best["tsd_present"],
        tsd_len=best["tsd_len"],
        tsd_seq=best["tsd_seq"],
        terminality_score=best["terminality_score"],
        complexity_score=best["complexity_score"],
        window=best["window"],
        confidence=confidence,
    )

    if logger:
        logger.info(
            "[TIR] %s len=%d identity=%.3f window=%d tsd=%s confidence=%s",
            family_id,
            det.tir_len,
            det.tir_identity,
            det.window,
            det.tsd_seq if det.tsd_present else "none",
            det.confidence,
        )

    return det


def _tir_worker(args):
    (
        family_id,
        seq,
        windows,
        min_len,
        max_len,
        min_identity,
        min_complexity,
    ) = args

    best_det = None

    for window in windows:
        det = best_terminal_inverted_repeat(
            seq=seq,
            family_id=family_id,
            window=window,
            min_len=min_len,
            max_len=max_len,
            min_identity=min_identity,
            min_complexity=min_complexity,
            logger=None,
        )

        if det is None:
            continue

        det = refine_tir_detection(det)

        if not det.tir_present:
            continue

        if (
            best_det is None
            or candidate_is_better(
                {
                    "tir_len": det.tir_len,
                    "tir_identity": det.tir_identity,
                    "terminality_score": det.terminality_score,
                    "tsd_present": det.tsd_present,
                    "window": det.window,
                },
                {
                    "tir_len": best_det.tir_len,
                    "tir_identity": best_det.tir_identity,
                    "terminality_score": best_det.terminality_score,
                    "tsd_present": best_det.tsd_present,
                    "window": best_det.window,
                },
            )
        ):
            best_det = det

    if best_det is None:
        return TirDetection(
            family_id=family_id,
            tir_present=False,
        )

    return best_det


def detect_tirs_from_fasta(
    fasta: str | Path,
    windows: list[int] | None = None,
    window: int | None = None,
    min_len: int = 12,
    max_len: int = 60,
    min_identity: float = 0.85,
    min_complexity: float = 0.50,
    threads: int = 1,
    logger=None,
) -> list[TirDetection]:
    if windows is None:
        if window is not None:
            windows = [window]
    
        else:
            #windows = [60, 90, 120, 180, 250]
            windows = [5, 10, 20, 60, 80]

    records = list(SeqIO.parse(str(fasta), "fasta"))
    total = len(records)

    jobs = [
        (
            clean_family_id(rec.id),
            str(rec.seq).upper(),
            windows,
            min_len,
            max_len,
            min_identity,
            min_complexity,
        )
        for rec in records
    ]

    detections: list[TirDetection] = []

    if logger:
        logger.info(
            "[TIR] starting parallel scan: families=%d threads=%d windows=%s",
            total,
            threads,
            ",".join(str(w) for w in windows),
        )

    if threads <= 1:
        for i, job in enumerate(jobs, start=1):
            det = _tir_worker(job)
            detections.append(det)

            if logger and (i == 1 or i % 25 == 0 or i == total):
                logger.info(
                    "[TIR] progress %d/%d families scanned; detected=%d",
                    i,
                    total,
                    sum(1 for d in detections if d.tir_present),
                )
    else:
        with ProcessPoolExecutor(max_workers=threads) as exe:
            for i, det in enumerate(exe.map(_tir_worker, jobs), start=1):
                detections.append(det)

                if logger and (i == 1 or i % 25 == 0 or i == total):
                    logger.info(
                        "[TIR] progress %d/%d families scanned; detected=%d",
                        i,
                        total,
                        sum(1 for d in detections if d.tir_present),
                    )

    detections.sort(key=lambda x: x.family_id)

    if logger:
        logger.info(
            "[TIR] final detected=%d total=%d",
            sum(1 for d in detections if d.tir_present),
            total,
        )

    return detections


def write_tir_structure_tsv(
    detections: list[TirDetection],
    outpath: str | Path,
) -> None:
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

        for det in detections:
            if not det.tir_present:
                continue

            source = (
                f"heuristic_tir:"
                f"len={det.tir_len};"
                f"window={det.window};"
                f"identity={det.tir_identity:.3f};"
                f"left={det.left_tir_start}-{det.left_tir_end};"
                f"right={det.right_tir_start}-{det.right_tir_end};"
                f"tsd_present={det.tsd_present};"
                f"tsd_len={det.tsd_len};"
                f"tsd_seq={det.tsd_seq};"
                f"terminality={det.terminality_score:.3f};"
                f"complexity={det.complexity_score:.3f};"
                f"confidence={det.confidence}"
            )

            writer.writerow({
                "family_id": det.family_id,
                "evidence_type": "TIR",
                "score": round(det.tir_identity, 3),
                "source": source,
            })

            if det.tsd_present:
                writer.writerow({
                    "family_id": det.family_id,
                    "evidence_type": "TSD",
                    "score": round(det.terminality_score, 3),
                    "source": (
                        f"heuristic_tsd:"
                        f"len={det.tsd_len};"
                        f"seq={det.tsd_seq};"
                        f"from_tir=true"
                    ),
                })

