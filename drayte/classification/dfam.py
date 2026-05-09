from __future__ import annotations

import os
import subprocess
import sys

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from typing import List

from .ids import clean_family_id
from drayte.utils.fasta_split import split_fasta_balanced_by_length
from concurrent.futures import ThreadPoolExecutor, as_completed

@dataclass
class DfamHit:
    family_id: str
    model_name: str
    accession: str
    evalue: float
    score: float
    ali_start: int
    ali_end: int
    dfam_class: str = "Unknown"
    dfam_order: str = "Unknown"
    dfam_superfamily: str = "Unknown"
    dfam_ct: str = ""

def run_nhmmer(
    dfam_db: str | Path,
    consensus_fasta: str | Path,
    tblout: str | Path,
    nhmmer_bin: str = "nhmmer",
    cpu: int = 1,
) -> Path:
    tblout = Path(tblout)
    tblout.parent.mkdir(parents=True, exist_ok=True)

    pid = os.getpid()
    start = datetime.now()

    cmd = [
        nhmmer_bin,
        "--cpu", str(cpu),
        "--tblout", str(tblout),
        str(dfam_db),
        str(consensus_fasta),
    ]

    with open(tblout.with_suffix(".log"), "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)

    end = datetime.now()

    return tblout

def parse_nhmmer_tblout(
    tblout: str | Path,
    max_evalue: float = 1e-5,
    min_score: float = 20.0,
    dfam_metadata: dict[str, dict[str, str]] | None = None,
) -> List[DfamHit]:

    hits: List[DfamHit] = []
    dfam_metadata = dfam_metadata or {}

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

            meta = dfam_metadata.get(model_name, {})

            hits.append(
                DfamHit(
                    family_id=family_id,
                    model_name=model_name,
                    accession=accession,
                    evalue=evalue,
                    score=score,
                    ali_start=ali_start,
                    ali_end=ali_end,
                    dfam_class=meta.get("dfam_class", "Unknown"),
                    dfam_order=meta.get("dfam_order", "Unknown"),
                    dfam_superfamily=meta.get(
                        "dfam_superfamily",
                        "Unknown",
                    ),
                    dfam_ct=meta.get("dfam_ct", ""),
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
    """
    Run nhmmer on a consensus FASTA split into length-balanced chunks.

    The merged output is written to:

        <outdir>/dfam.merged.tblout

    Chunk FASTAs are written to:

        <outdir>/chunks/

    Per-chunk tblout files are written to:

        <outdir>/tblouts/
    """

    dfam_db = Path(dfam_db)
    consensus_fasta = Path(consensus_fasta)
    outdir = Path(outdir)

    chunks_dir = outdir / "chunks"
    tblout_dir = outdir / "tblouts"

    chunks_dir.mkdir(parents=True, exist_ok=True)
    tblout_dir.mkdir(parents=True, exist_ok=True)

    if chunks < 1:
        chunks = 1

    fasta_chunks = split_fasta_balanced_by_length(
        input_fasta=consensus_fasta,
        outdir=chunks_dir,
        n_chunks=chunks,
        prefix="chunk",
    )

    if not fasta_chunks:
        merged = outdir / "dfam.merged.tblout"
        merged.write_text("")
        return merged

    if max_parallel is None:
        max_parallel = len(fasta_chunks)

    max_parallel = max(1, min(max_parallel, len(fasta_chunks)))

    def _run_one(chunk_fasta: Path) -> Path:
        tblout = tblout_dir / f"{chunk_fasta.stem}.tblout"

        return run_nhmmer(
            dfam_db=dfam_db,
            consensus_fasta=chunk_fasta,
            tblout=tblout,
            nhmmer_bin=nhmmer_bin,
            cpu=cpu_per_job,
        )

    tblouts: list[Path] = []

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        futures = {
            executor.submit(_run_one, chunk): chunk
            for chunk in fasta_chunks
        }

        total = len(futures)
        completed = 0

        for future in as_completed(futures):
            chunk = futures[future]
            tblout = future.result()
            tblouts.append(tblout)

            completed += 1

            print(
                f"[DRayTE][Dfam] progress "
                f"{completed}/{total} chunks completed "
                f"| last={chunk.name}",
                file=sys.stderr,
                flush=True,
            )

    merged = outdir / "dfam.merged.tblout"

    with open(merged, "w") as out:
        for tblout in sorted(tblouts):
            with open(tblout) as fh:
                for line in fh:
                    out.write(line)

    return merged
