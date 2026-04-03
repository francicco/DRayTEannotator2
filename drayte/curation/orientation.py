from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, Set

from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement_fasta(infile: Path, outfile: Path) -> None:
    records = []
    for rec in SeqIO.parse(str(infile), "fasta"):
        rec.seq = rec.seq.reverse_complement()
        records.append(rec)
    with open(outfile, "w") as handle:
        SeqIO.write(records, handle, "fasta")


def run_blastx_orientation(
    family_fa: Path,
    db_dmnd: Path,
    outfile: Path,
    diamond_bin: str,
    threads: int,
    logger,
) -> Path:
    if outfile.exists() and outfile.stat().st_size > 0:
        return outfile

    cmd = [
        diamond_bin,
        "blastx",
        "--threads", str(threads),
        "--db", str(db_dmnd),
        "--query", str(family_fa),
        "--outfmt", "6",
        "--evalue", "1e-15",
        "--out", str(outfile),
    ]
    logger.info("Running blastx orientation check for %s", family_fa.name)
    subprocess.run(cmd, check=True)
    return outfile


def maybe_reorient_family(
    short_id: str,
    family_dir: Path,
    db_dmnd: Path,
    diamond_bin: str,
    threads: int,
    logger,
) -> bool:
    rep = family_dir / f"{short_id}_rep.fa"
    rep_mod = family_dir / f"{short_id}_rep_mod.fa"
    msa = family_dir / f"{short_id}_MSA_extended.fa"
    out = family_dir / f"{short_id}_extended_rep_blastx.out"

    if not rep.exists():
        return False

    run_blastx_orientation(rep, db_dmnd, out, diamond_bin, threads, logger)

    if not out.exists() or out.stat().st_size == 0:
        return False

    with open(out) as handle:
        first = handle.readline().strip().split("\t")
    if len(first) < 8:
        return False

    start = int(first[6])
    end = int(first[7])

    if start <= end:
        return False

    logger.info("Reverse-complementing %s", short_id)

    tmp = rep.with_suffix(".tmp")
    reverse_complement_fasta(rep, tmp)
    tmp.replace(rep)

    if msa.exists():
        tmp = msa.with_suffix(".tmp")
        reverse_complement_fasta(msa, tmp)
        tmp.replace(msa)

    if rep_mod.exists():
        tmp = rep_mod.with_suffix(".tmp")
        reverse_complement_fasta(rep_mod, tmp)
        tmp.replace(rep_mod)

    return True
