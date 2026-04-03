from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Dict, Optional


def ensure_repeatpeps_db(
    db_dir: Path,
    diamond_bin: str,
    logger,
) -> Path:
    db_dir.mkdir(parents=True, exist_ok=True)

    repeatpeps_lib = db_dir / "RepeatPeps.lib"
    repeatpeps_dmnd = db_dir / "RepeatPeps.lib.dmnd"

    if not repeatpeps_lib.exists():
        raise FileNotFoundError(
            f"RepeatPeps.lib not found: {repeatpeps_lib}. "
            "Download or place it there before running."
        )

    if not repeatpeps_dmnd.exists():
        logger.info("Building DIAMOND database for RepeatPeps")
        subprocess.run(
            [diamond_bin, "makedb", "--in", str(repeatpeps_lib), "--db", str(repeatpeps_lib)],
            check=True,
        )

    return repeatpeps_dmnd


def run_diamond_blastp(
    query_fasta: Path,
    db_path: Path,
    out_tsv: Path,
    diamond_bin: str,
    threads: int,
    logger,
) -> Path:
    if out_tsv.exists() and out_tsv.stat().st_size > 0:
        logger.info("DIAMOND blastp already exists: %s", out_tsv)
        return out_tsv

    logger.info("Running DIAMOND blastp against RepeatPeps")
    subprocess.run(
        [
            diamond_bin,
            "blastp",
            "--threads", str(threads),
            "--db", str(db_path),
            "--query", str(query_fasta),
            "--out", str(out_tsv),
            "--outfmt", "6",
            "--evalue", "1e-15",
        ],
        check=True,
    )
    return out_tsv


def run_diamond_blastx(
    query_fasta: Path,
    db_path: Path,
    out_tsv: Path,
    diamond_bin: str,
    threads: int,
    logger,
) -> Path:
    if out_tsv.exists() and out_tsv.stat().st_size > 0:
        logger.info("DIAMOND blastx already exists: %s", out_tsv)
        return out_tsv

    logger.info("Running DIAMOND blastx against RepeatPeps: %s", query_fasta.name)
    subprocess.run(
        [
            diamond_bin,
            "blastx",
            "--threads", str(threads),
            "--db", str(db_path),
            "--query", str(query_fasta),
            "--out", str(out_tsv),
            "--outfmt", "6",
            "--evalue", "1e-15",
        ],
        check=True,
    )
    return out_tsv


def parse_top_hit_per_query(blast_tsv: Path) -> Dict[str, dict]:
    """
    Return best hit per query by highest bitscore, then longest alignment.
    DIAMOND outfmt 6 default columns assumed:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    """
    best: Dict[str, dict] = {}

    with open(blast_tsv) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            qseqid = row[0]
            hit = {
                "qseqid": row[0],
                "sseqid": row[1],
                "pident": float(row[2]),
                "length": int(row[3]),
                "qstart": int(row[6]),
                "qend": int(row[7]),
                "sstart": int(row[8]),
                "send": int(row[9]),
                "evalue": float(row[10]),
                "bitscore": float(row[11]),
            }

            if qseqid not in best:
                best[qseqid] = hit
            else:
                old = best[qseqid]
                if (hit["bitscore"], hit["length"]) > (old["bitscore"], old["length"]):
                    best[qseqid] = hit

    return best
