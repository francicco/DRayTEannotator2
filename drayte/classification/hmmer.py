from __future__ import annotations

import subprocess
from dataclasses import dataclass

from .ids import clean_family_id
from pathlib import Path
from typing import Dict, List

from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO

@dataclass
class DomainHit:
    family_id: str
    orf_id: str
    domain: str
    accession: str
    evalue: float
    score: float
    hmm_start: int
    hmm_end: int
    ali_start: int
    ali_end: int

def split_fasta(
    fasta: str | Path,
    outdir: str | Path,
    chunks: int,
    suffix: str = ".fa",
) -> list[Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(str(fasta), "fasta"))

    if not records:
        return []

    if chunks <= 1:
        out = outdir / f"chunk_001{suffix}"
        SeqIO.write(records, str(out), "fasta")
        return [out]

    chunk_size = max(1, (len(records) + chunks - 1) // chunks)
    chunk_files = []

    for i in range(0, len(records), chunk_size):
        chunk_records = records[i:i + chunk_size]
        chunk_id = len(chunk_files) + 1
        out = outdir / f"chunk_{chunk_id:03d}{suffix}"
        SeqIO.write(chunk_records, str(out), "fasta")
        chunk_files.append(out)

    return chunk_files

def merge_domtblouts(
    domtblouts: list[Path],
    merged: str | Path,
) -> Path:
    merged = Path(merged)
    merged.parent.mkdir(parents=True, exist_ok=True)

    with open(merged, "w") as out:
        wrote_header = False

        for path in domtblouts:
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        if not wrote_header:
                            out.write(line)
                        continue

                    out.write(line)

            wrote_header = True

    return merged

def run_hmmscan_parallel_chunks(
    proteins_fasta: str | Path,
    hmm_db: str | Path,
    outdir: str | Path,
    chunks: int,
    hmmscan_bin: str = "hmmscan",
    cpu_per_job: int = 1,
    max_parallel: int | None = None,
) -> Path:
    outdir = Path(outdir)

    chunk_dir = outdir / "chunks"
    scan_dir = outdir / "hmmscan_chunks"
    scan_dir.mkdir(parents=True, exist_ok=True)

    chunk_files = split_fasta(
        proteins_fasta,
        chunk_dir,
        chunks,
        suffix=".fa",
    )

    if not chunk_files:
        merged = outdir / "domains.domtblout"
        merged.write_text("")
        return merged

    if max_parallel is None:
        max_parallel = len(chunk_files)

    max_parallel = max(1, min(max_parallel, len(chunk_files)))

    def _run_one(chunk_fasta: Path) -> Path:
        domtblout = scan_dir / f"{chunk_fasta.stem}.domtblout"

        return run_hmmscan(
            hmm_db=hmm_db,
            proteins_fasta=chunk_fasta,
            domtblout=domtblout,
            hmmscan_bin=hmmscan_bin,
            cpu=cpu_per_job,
        )

    domtblouts = []

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        futures = {
            executor.submit(_run_one, chunk): chunk
            for chunk in chunk_files
        }

        for future in as_completed(futures):
            domtblouts.append(future.result())

    merged = outdir / "domains.domtblout"

    return merge_domtblouts(
        sorted(domtblouts),
        merged,
    )

def run_hmmscan(
    hmm_db: str | Path,
    proteins_fasta: str | Path,
    domtblout: str | Path,
    hmmscan_bin: str = "hmmscan",
    cpu: int = 1,
) -> Path:
    domtblout = Path(domtblout)
    domtblout.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        hmmscan_bin,
        "--cpu", str(cpu),
        "--domtblout", str(domtblout),
        str(hmm_db),
        str(proteins_fasta),
    ]

    with open(domtblout.with_suffix(".log"), "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)

    return domtblout


def family_from_orf_id(orf_id: str) -> str:
    if "_orf" in orf_id:
        return clean_family_id(orf_id.rsplit("_orf", 1)[0])
    return clean_family_id(orf_id)


def parse_domtblout(
    domtblout: str | Path,
    max_evalue: float = 1e-5,
    min_score: float = 20.0,
) -> List[DomainHit]:
    hits: List[DomainHit] = []

    with open(domtblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split()
            if len(parts) < 23:
                continue

            domain = parts[0]
            accession = parts[1]
            orf_id = parts[3]

            evalue = float(parts[12])
            score = float(parts[13])

            if evalue > max_evalue:
                continue
            if score < min_score:
                continue

            hits.append(
                DomainHit(
                    family_id=family_from_orf_id(orf_id),
                    orf_id=orf_id,
                    domain=domain,
                    accession=accession,
                    evalue=evalue,
                    score=score,
                    hmm_start=int(parts[15]),
                    hmm_end=int(parts[16]),
                    ali_start=int(parts[17]),
                    ali_end=int(parts[18]),
                )
            )

    return hits


def summarize_domains_by_family(hits: List[DomainHit]) -> Dict[str, dict]:
    summary: Dict[str, dict] = {}

    for hit in hits:
        fam = summary.setdefault(
            hit.family_id,
            {
                "domains": set(),
                "best_domain_score": 0.0,
                "domain_hits": 0,
            },
        )

        fam["domains"].add(hit.domain)
        fam["best_domain_score"] = max(fam["best_domain_score"], hit.score)
        fam["domain_hits"] += 1

    for fam in summary.values():
        fam["domains"] = sorted(fam["domains"])

    return summary
