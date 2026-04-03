
from __future__ import annotations

import subprocess
from pathlib import Path


def run_orf_detection(input_fasta: Path, outdir: Path, getorf_bin: str, logger) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)

    orf_fasta = outdir / f"{input_fasta.stem}.orfs.fa"

    if not orf_fasta.exists() or orf_fasta.stat().st_size == 0:
        cmd = [
            getorf_bin,
            "-sequence", str(input_fasta),
            "-outseq", str(orf_fasta),
            "-minsize", "300",
        ]
        logger.info("Running getorf")
        subprocess.run(cmd, check=True)
    else:
        logger.info("ORF output already exists: %s", orf_fasta)

    return {"orf_fasta": str(orf_fasta)}
