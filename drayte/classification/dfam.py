from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List

from .ids import clean_family_id


@dataclass
class DfamHit:
    family_id: str
    model_name: str
    accession: str
    evalue: float
    score: float
    ali_start: int
    ali_end: int


def run_nhmmer(
    dfam_db: str | Path,
    consensus_fasta: str | Path,
    tblout: str | Path,
    nhmmer_bin: str = "nhmmer",
    cpu: int = 1,
) -> Path:
    tblout = Path(tblout)
    tblout.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        nhmmer_bin,
        "--cpu", str(cpu),
        "--tblout", str(tblout),
        str(dfam_db),
        str(consensus_fasta),
    ]

    with open(tblout.with_suffix(".log"), "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)

    return tblout


def parse_nhmmer_tblout(
    tblout: str | Path,
    max_evalue: float = 1e-5,
    min_score: float = 20.0,
) -> List[DfamHit]:
    hits: List[DfamHit] = []

    with open(tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split()
            if len(parts) < 14:
                continue

            family_id = clean_family_id(parts[0])
            model_name = parts[2]
            accession = parts[3]
            ali_start = int(parts[6])
            ali_end = int(parts[7])
            evalue = float(parts[12])
            score = float(parts[13])

            if evalue > max_evalue:
                continue

            if score < min_score:
                continue

            hits.append(
                DfamHit(
                    family_id=family_id,
                    model_name=model_name,
                    accession=accession,
                    evalue=evalue,
                    score=score,
                    ali_start=ali_start,
                    ali_end=ali_end,
                )
            )

    return hits

def best_dfam_hits_by_family(hits: List[DfamHit]) -> dict[str, DfamHit]:
    best = {}

    for hit in hits:
        current = best.get(hit.family_id)

        if current is None:
            best[hit.family_id] = hit
            continue

        if hit.score > current.score:
            best[hit.family_id] = hit
            continue

        if hit.score == current.score and hit.evalue < current.evalue:
            best[hit.family_id] = hit

    return best


def split_fasta_round_robin(
    fasta: str | Path,
    outdir: str | Path,
    chunks: int,
) -> list[Path]:
    from Bio import SeqIO

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    handles = []
    paths = []

    for i in range(chunks):
        path = outdir / f"chunk_{i + 1:03d}.fa"
        paths.append(path)
        handles.append(open(path, "w"))

    try:
        for idx, rec in enumerate(SeqIO.parse(str(fasta), "fasta")):
            handle = handles[idx % chunks]
            SeqIO.write(rec, handle, "fasta")
    finally:
        for handle in handles:
            handle.close()

    return [
        p for p in paths
        if p.exists() and p.stat().st_size > 0
    ]


def merge_tblout_files(
    tblouts: list[Path],
    merged_tblout: str | Path,
) -> Path:
    merged_tblout = Path(merged_tblout)
    merged_tblout.parent.mkdir(parents=True, exist_ok=True)

    with open(merged_tblout, "w") as out:
        for path in tblouts:
            with open(path) as fh:
                for line in fh:
                    out.write(line)

    return merged_tblout


def run_nhmmer_parallel_chunks(
    dfam_db: str | Path,
    consensus_fasta: str | Path,
    outdir: str | Path,
    chunks: int,
    nhmmer_bin: str = "nhmmer",
    cpu_per_job: int = 1,
    max_parallel: int | None = None,
) -> Path:
    from concurrent.futures import ThreadPoolExecutor, as_completed

    outdir = Path(outdir)
    chunks_dir = outdir / "chunks"
    tblout_dir = outdir / "tblouts"

    chunks_dir.mkdir(parents=True, exist_ok=True)
    tblout_dir.mkdir(parents=True, exist_ok=True)

    fasta_chunks = split_fasta_round_robin(
        fasta=consensus_fasta,
        outdir=chunks_dir,
        chunks=chunks,
    )

    if max_parallel is None:
        max_parallel = chunks

    tblouts = []

    def _run_one(chunk_fa: Path) -> Path:
        tblout = tblout_dir / f"{chunk_fa.stem}.tblout"

        run_nhmmer(
            dfam_db=dfam_db,
            consensus_fasta=chunk_fa,
            tblout=tblout,
            nhmmer_bin=nhmmer_bin,
            cpu=cpu_per_job,
        )

        return tblout

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        futures = [
            executor.submit(_run_one, chunk)
            for chunk in fasta_chunks
        ]

        for future in as_completed(futures):
            tblouts.append(future.result())

    tblouts = sorted(tblouts)

    merged = outdir / "dfam.merged.tblout"

    return merge_tblout_files(
        tblouts=tblouts,
        merged_tblout=merged,
    )
