from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO

from .repeatpeps import ensure_repeatpeps_db, run_diamond_blastp, parse_top_hit_per_query


def normalize_rc_header(header: str) -> tuple[str, str, str]:
    """
    Convert RepeatClassifier-style header into:
    original_header, name, class, family
    """
    original = header.strip()

    if "#" in original:
        name, clf = original.split("#", 1)
    else:
        name, clf = original, "Unknown/Unknown"

    clf = clf.strip()

    replacements = {
        "Unknown": "Unknown/Unknown",
        "Satellite": "Satellite/Satellite",
        "LTR ": "LTR/Unknown",
        "DNA ": "DNA/Unknown",
        "tRNA ": "tRNA/Nothing",
        "LINE ": "LINE/Unknown",
    }
    clf = replacements.get(clf, clf)

    if "/" in clf:
        te_class, family = clf.split("/", 1)
    else:
        te_class, family = clf, "Unknown"

    return original, name, te_class, family


def run_getorf(
    input_fasta: Path,
    output_fasta: Path,
    getorf_bin: str,
    min_orf: int,
    logger,
) -> Path:
    if output_fasta.exists() and output_fasta.stat().st_size > 0:
        logger.info("getorf output already exists: %s", output_fasta)
        return output_fasta

    logger.info("Running getorf on %s", input_fasta.name)
    subprocess.run(
        [
            getorf_bin,
            "-sequence", str(input_fasta),
            "-outseq", str(output_fasta),
            "-minsize", str(min_orf),
        ],
        check=True,
    )
    return output_fasta


def summarize_orf_hits_by_family(
    blastp_tsv: Path,
) -> Dict[str, dict]:
    """
    Group ORF blastp hits back to consensus family using the ORF header prefix.
    Assumes getorf-derived IDs begin with the original sequence name.
    """
    best_orf_hits = parse_top_hit_per_query(blastp_tsv)
    by_family: Dict[str, List[dict]] = {}

    for qid, hit in best_orf_hits.items():
        family = qid.split("_", 1)[0]
        by_family.setdefault(family, []).append(hit)

    summary: Dict[str, dict] = {}
    for family, hits in by_family.items():
        hits = sorted(hits, key=lambda x: (x["bitscore"], x["length"]), reverse=True)
        top = hits[0]
        subject = top["sseqid"].replace("--", "#")
        if "#" in subject:
            _, rest = subject.split("#", 1)
        else:
            rest = subject

        if "/" in rest:
            hit_class, hit_family = rest.split("/", 1)
        else:
            hit_class, hit_family = rest, "NOHIT"

        summary[family] = {
            "orf_hit_count": len(hits),
            "top_orf_class": hit_class,
            "top_orf_family": hit_family,
            "top_orf_subject": subject,
            "top_orf_align_len": top["length"],
            "top_orf_bitscore": top["bitscore"],
        }

    return summary


def build_family_table_from_classified_library(
    classified_library: Path,
    outdir: Path,
    getorf_bin: str,
    diamond_bin: str,
    repeatpeps_db_dir: Path,
    min_orf: int,
    threads: int,
    logger,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)

    orf_fasta = outdir / f"{classified_library.name}_getorf.fa"
    run_getorf(classified_library, orf_fasta, getorf_bin, min_orf, logger)

    repeatpeps_db = ensure_repeatpeps_db(repeatpeps_db_dir, diamond_bin, logger)
    blastp_tsv = outdir / f"{classified_library.name}_rep_blastp.out"
    run_diamond_blastp(orf_fasta, repeatpeps_db, blastp_tsv, diamond_bin, threads, logger)

    orf_summary = summarize_orf_hits_by_family(blastp_tsv)

    rows = []
    for rec in SeqIO.parse(str(classified_library), "fasta"):
        original, name, te_class, family = normalize_rc_header(rec.id)
        hit = orf_summary.get(name, None)

        row = {
            "original_header": original,
            "name": name,
            "class": te_class,
            "family": family,
            "consensus_len": len(rec.seq),
            "orf_hit_count": hit["orf_hit_count"] if hit else 0,
            "top_orf_class": hit["top_orf_class"] if hit else "NOHIT",
            "top_orf_family": hit["top_orf_family"] if hit else "NOHIT",
            "top_orf_subject": hit["top_orf_subject"] if hit else "NOHIT",
            "top_orf_align_len": hit["top_orf_align_len"] if hit else 0,
            "top_orf_bitscore": hit["top_orf_bitscore"] if hit else 0.0,
        }
        rows.append(row)

    table = outdir / "family_table.tsv"
    with open(table, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "original_header",
                "name",
                "class",
                "family",
                "consensus_len",
                "orf_hit_count",
                "top_orf_class",
                "top_orf_family",
                "top_orf_subject",
                "top_orf_align_len",
                "top_orf_bitscore",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Wrote family table: %s", table)
    return table
