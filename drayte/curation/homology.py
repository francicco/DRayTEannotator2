from __future__ import annotations

import subprocess
from pathlib import Path


def run_repeatpep_search(
    input_fasta: Path,
    outdir: Path,
    diamond_bin: str,
    repeatpeps_db: str | None,
    threads: int,
    logger,
) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)

    if not repeatpeps_db:
        logger.info("No RepeatPeps database configured; skipping homology search")
        return {"diamond_tsv": None}

    diamond_tsv = outdir / f"{input_fasta.stem}.diamond.tsv"

    if not diamond_tsv.exists() or diamond_tsv.stat().st_size == 0:
        cmd = [
            diamond_bin,
            "blastx",
            "--db", repeatpeps_db,
            "--query", str(input_fasta),
            "--out", str(diamond_tsv),
            "--outfmt", "6",
            "--threads", str(threads),
        ]
        logger.info("Running DIAMOND blastx against RepeatPeps")
        subprocess.run(cmd, check=True)
    else:
        logger.info("Homology output already exists: %s", diamond_tsv)

    return {"diamond_tsv": str(diamond_tsv)}
