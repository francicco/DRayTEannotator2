from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO


def run_getorf(input_fasta: str | Path, outdir: str | Path, getorf_bin: str, min_orf: int, logger) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_fa = outdir / f"{Path(input_fasta).name}_getorf.fa"

    if out_fa.exists() and out_fa.stat().st_size > 0:
        logger.info("getorf output already exists: %s", out_fa)
        return out_fa

    cmd = [
        getorf_bin,
        str(input_fasta),
        str(out_fa),
        "-minsize",
        str(min_orf),
    ]
    logger.info("Running getorf")
    subprocess.run(cmd, check=True)
    return out_fa


def collect_orf_lengths(getorf_fasta: str | Path) -> Dict[str, List[int]]:
    result: Dict[str, List[int]] = {}
    for record in SeqIO.parse(str(getorf_fasta), "fasta"):
        parent = record.id.split("_")[0]
        result.setdefault(parent, []).append(len(record.seq))
    for k in result:
        result[k] = sorted(result[k], reverse=True)
    return result
