from __future__ import annotations

import csv
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass
class FamilyDecision:
    family_id: str
    te_class: str
    te_family: str
    consensus_len: int
    orf_hit_count: int
    top_orf_class: str
    top_orf_family: str
    top_orf_subject: str
    top_orf_align_len: int
    top_orf_bitscore: float
    score: float
    decision: str
    reason: str

def normalize_repeatmasker_classification(classification: str) -> str:
    classification = classification.strip()

    # Penelope must be its own top-level class for proper RepeatMasker summary
    if classification.startswith("LINE/Penelope"):
        return "LINE/Penelope"

    # Helitrons should be reported as rolling-circle elements
    if classification in {"DNA/Helitron", "Helitron"}:
        return "RC/Helitron"

    if classification.startswith("RC/Helitron"):
        return "RC/Helitron"

    return classification

def load_family_table(family_table_tsv: Path) -> List[dict]:
    with open(family_table_tsv) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def score_family(row: dict) -> tuple[float, str, str]:
    """
    Conservative first-pass scoring.
    Output:
      score, decision, reason
    """

    consensus_len = int(row["consensus_len"])
    orf_hit_count = int(row["orf_hit_count"])
    top_orf_bitscore = float(row["top_orf_bitscore"])
    top_orf_class = row["top_orf_class"]
    te_class = row["class"]

    score = 0.0
    reasons: List[str] = []

    # length
    if consensus_len >= 4000:
        score += 3
        reasons.append("long_consensus")
    elif consensus_len >= 1000:
        score += 2
        reasons.append("medium_consensus")
    elif consensus_len >= 300:
        score += 1
        reasons.append("short_consensus")
    else:
        score -= 2
        reasons.append("very_short")

    # ORF evidence
    if orf_hit_count >= 3:
        score += 3
        reasons.append("multiple_orf_hits")
    elif orf_hit_count >= 1:
        score += 1.5
        reasons.append("some_orf_hits")
    else:
        reasons.append("no_orf_hits")

    # protein similarity strength
    if top_orf_bitscore >= 300:
        score += 3
        reasons.append("strong_protein_hit")
    elif top_orf_bitscore >= 100:
        score += 2
        reasons.append("moderate_protein_hit")
    elif top_orf_bitscore > 0:
        score += 1
        reasons.append("weak_protein_hit")

    # consistency between RepeatClassifier and protein evidence
    if top_orf_class != "NOHIT":
        if te_class == top_orf_class:
            score += 2
            reasons.append("class_consistent")
        elif te_class == "Unknown":
            score += 1
            reasons.append("protein_resolves_unknown")
        else:
            score -= 1
            reasons.append("class_conflict")

    # hard flags
    if consensus_len < 200 and orf_hit_count == 0:
        return score, "discard", "too_short_no_support"

    if score >= 6:
        return score, "keep", ",".join(reasons)
    if score >= 3:
        return score, "review", ",".join(reasons)
    return score, "discard", ",".join(reasons)


def make_decisions(rows: List[dict]) -> List[FamilyDecision]:
    decisions: List[FamilyDecision] = []

    for row in rows:
        score, decision, reason = score_family(row)
        decisions.append(
            FamilyDecision(
                family_id=row["name"],
                te_class=row["class"],
                te_family=row["family"],
                consensus_len=int(row["consensus_len"]),
                orf_hit_count=int(row["orf_hit_count"]),
                top_orf_class=row["top_orf_class"],
                top_orf_family=row["top_orf_family"],
                top_orf_subject=row["top_orf_subject"],
                top_orf_align_len=int(row["top_orf_align_len"]),
                top_orf_bitscore=float(row["top_orf_bitscore"]),
                score=score,
                decision=decision,
                reason=reason,
            )
        )

    return decisions


def load_fasta_as_dict(fasta: Path) -> Dict[str, SeqRecord]:
    return {rec.id: rec for rec in SeqIO.parse(str(fasta), "fasta")}


def normalize_record_header_for_repeatmasker(rec: SeqRecord) -> SeqRecord:
    rec = rec[:]

    if "#" in rec.id:
        name, classification = rec.id.split("#", 1)
        classification = normalize_repeatmasker_classification(classification)
        rec.id = f"{name}#{classification}"
        rec.name = rec.id
        rec.description = rec.id

    return rec

def write_curated_full_library(
    classified_library: Path,
    decisions: List[FamilyDecision],
    outfile: Path,
) -> int:
    records = load_fasta_as_dict(classified_library)
    keep_ids = {d.family_id for d in decisions if d.decision in {"keep", "review"}}

    selected = [
        normalize_record_header_for_repeatmasker(records[k])
        for k in keep_ids
        if k in records
    ]

    with open(outfile, "w") as handle:
        SeqIO.write(selected, handle, "fasta")

    return len(selected)

def write_metadata(decisions: List[FamilyDecision], outfile: Path) -> None:
    with open(outfile, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "family_id",
            "class",
            "family",
            "consensus_len",
            "orf_hit_count",
            "top_orf_class",
            "top_orf_family",
            "top_orf_subject",
            "top_orf_align_len",
            "top_orf_bitscore",
            "score",
            "decision",
            "reason",
        ])
        for d in decisions:
            writer.writerow([
                d.family_id,
                d.te_class,
                d.te_family,
                d.consensus_len,
                d.orf_hit_count,
                d.top_orf_class,
                d.top_orf_family,
                d.top_orf_subject,
                d.top_orf_align_len,
                d.top_orf_bitscore,
                d.score,
                d.decision,
                d.reason,
            ])


def write_run_summary(decisions: List[FamilyDecision], outfile: Path, n_full: int, n_nr: int) -> None:
    keep = sum(1 for d in decisions if d.decision == "keep")
    review = sum(1 for d in decisions if d.decision == "review")
    discard = sum(1 for d in decisions if d.decision == "discard")

    with open(outfile, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["n_total", len(decisions)])
        writer.writerow(["n_keep", keep])
        writer.writerow(["n_review", review])
        writer.writerow(["n_discard", discard])
        writer.writerow(["n_curated_full", n_full])
        writer.writerow(["n_curated_nr", n_nr])


def run_cd_hit_est(
    input_fasta: Path,
    output_fasta: Path,
    cdhit_bin: str,
    identity: float,
    threads: int,
    logger,
) -> Path:
    clstr = Path(str(output_fasta) + ".clstr")

    if output_fasta.exists() and output_fasta.stat().st_size > 0:
        logger.info("Non-redundant library already exists: %s", output_fasta)
        return output_fasta

    cmd = [
        cdhit_bin,
        "-i", str(input_fasta),
        "-o", str(output_fasta),
        "-c", str(identity),
        "-G", "0",
        "-aS", "0.8",
        "-aL", "0.8",
        "-g", "1",
        "-M", "0",
        "-T", str(threads),
        "-d", "0",
    ]

    logger.info("Running CD-HIT-EST for non-redundant library")
    subprocess.run(cmd, check=True)

    return output_fasta


def finalize_curated_outputs(
    classified_library: Path,
    family_table_tsv: Path,
    outdir: Path,
    cdhit_bin: str | None,
    cdhit_identity: float,
    threads: int,
    logger,
) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)

    rows = load_family_table(family_table_tsv)
    decisions = make_decisions(rows)

    curated_full = outdir / "curated_full.fa"
    curated_metadata = outdir / "curated_metadata.tsv"
    curated_nr = outdir / "curated_nr.fa"
    run_summary = outdir / "run_summary.tsv"

    n_full = write_curated_full_library(classified_library, decisions, curated_full)
    write_metadata(decisions, curated_metadata)

    if cdhit_bin:
        run_cd_hit_est(
            input_fasta=curated_full,
            output_fasta=curated_nr,
            cdhit_bin=cdhit_bin,
            identity=cdhit_identity,
            threads=threads,
            logger=logger,
        )
        n_nr = sum(1 for _ in SeqIO.parse(str(curated_nr), "fasta"))
    else:
        logger.warning("No CD-HIT-EST configured; using curated_full.fa as curated_nr.fa")
        if not curated_nr.exists():
            curated_nr.write_text(curated_full.read_text())
        n_nr = n_full

    write_run_summary(decisions, run_summary, n_full=n_full, n_nr=n_nr)

    return {
        "curated_full": str(curated_full),
        "curated_nr": str(curated_nr),
        "curated_metadata": str(curated_metadata),
        "run_summary": str(run_summary),
    }
