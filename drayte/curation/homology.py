from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple


def ensure_repeatpeps_db(db_dir: str | Path, diamond_bin: str, logger) -> Path:
    db_dir = Path(db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)

    lib = db_dir / "RepeatPeps.lib"
    dmnd = db_dir / "RepeatPeps.lib.dmnd"

    if not lib.exists():
        logger.info("Downloading RepeatPeps.lib")
        subprocess.run(
            [
                "wget",
                "-O",
                str(lib),
                "https://raw.githubusercontent.com/rmhubley/RepeatMasker/master/Libraries/RepeatPeps.lib",
            ],
            check=True,
        )

    if not dmnd.exists():
        logger.info("Building DIAMOND database for RepeatPeps")
        subprocess.run(
            [diamond_bin, "makedb", "--in", str(lib), "--db", str(db_dir / "RepeatPeps.lib")],
            check=True,
        )

    return dmnd


def run_diamond_blastp(
    query_fasta: str | Path,
    db_dmnd: str | Path,
    outfile: str | Path,
    diamond_bin: str,
    threads: int,
    logger,
) -> Path:
    outfile = Path(outfile)
    if outfile.exists() and outfile.stat().st_size > 0:
        logger.info("DIAMOND blastp output already exists: %s", outfile)
        return outfile

    cmd = [
        diamond_bin,
        "blastp",
        "--threads", str(threads),
        "--db", str(db_dmnd),
        "--query", str(query_fasta),
        "--outfmt", "6",
        "--evalue", "1e-15",
        "--out", str(outfile),
    ]
    logger.info("Running DIAMOND blastp")
    subprocess.run(cmd, check=True)
    return outfile


def best_hits_per_orf(blastp_tsv: str | Path) -> Dict[str, Tuple[str, int]]:
    best: Dict[str, Tuple[str, int, float]] = {}
    with open(blastp_tsv) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            qseqid = row[0]
            sseqid = row[1].replace("--", "#")
            alen = int(row[3])
            bitscore = float(row[11])
            if qseqid not in best or bitscore > best[qseqid][2]:
                best[qseqid] = (sseqid, alen, bitscore)
    return {k: (v[0], v[1]) for k, v in best.items()}
