from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List

from Bio import SeqIO

from .homology import homology_from_repeatmasker_header
from .ids import clean_family_id


@dataclass
class MMseqsHit:
    query_id: str
    target_id: str
    target_header: str
    identity: float
    aln_len: int
    query_cov: float
    target_cov: float
    qstart: int
    qend: int
    tstart: int
    tend: int
    evalue: float
    bits: float


@dataclass
class RescueEvidence:
    family_id: str
    rescue_class: str
    rescue_order: str
    rescue_superfamily: str
    rescue_identity: float
    rescue_aln_len: int
    rescue_query_cov: float
    rescue_target_cov: float
    rescue_bits: float
    rescue_target: str


def sequence_lengths_by_clean_id(
    fasta: str | Path,
) -> dict[str, int]:
    return {
        clean_family_id(rec.id): len(rec.seq)
        for rec in SeqIO.parse(str(fasta), "fasta")
    }


def run_mmseqs_rescue(
    query_fasta: str | Path,
    reference_fasta: str | Path,
    outdir: str | Path,
    mmseqs_bin: str = "mmseqs",
    threads: int = 1,
    min_seq_id: float = 0.75,
    coverage: float = 0.50,
) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out_tsv = outdir / "mmseqs_rescue.tsv"
    tmpdir = outdir / "tmp"

    cmd = [
        mmseqs_bin,
        "easy-search",
        str(query_fasta),
        str(reference_fasta),
        str(out_tsv),
        str(tmpdir),
        "--search-type",
        "3",
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--cov-mode",
        "0",
        "--threads",
        str(threads),
    ]

    with open(outdir / "mmseqs_rescue.log", "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)

    return out_tsv


def parse_mmseqs_tsv(
    path: str | Path,
    query_lengths: dict[str, int] | None = None,
    target_lengths: dict[str, int] | None = None,
    min_identity: float = 0.80,
    min_aln_len: int = 500,
    min_query_cov: float = 0.0,
    min_target_cov: float = 0.0,
) -> List[MMseqsHit]:

    hits: List[MMseqsHit] = []

    query_lengths = query_lengths or {}
    target_lengths = target_lengths or {}

    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue

            parts = line.rstrip().split("\t")

            if len(parts) < 12:
                continue

            query_raw = parts[0]
            target_raw = parts[1]

            query_id = clean_family_id(query_raw)
            target_id = clean_family_id(target_raw)

            if query_id == target_id:
                continue

            identity = float(parts[2])
            aln_len = int(parts[3])

            qlen = query_lengths.get(query_id)
            tlen = target_lengths.get(target_id)

            query_cov = aln_len / qlen if qlen else 0.0
            target_cov = aln_len / tlen if tlen else 0.0

            if identity < min_identity:
                continue

            if aln_len < min_aln_len:
                continue

            if qlen and query_cov < min_query_cov:
                continue

            if tlen and target_cov < min_target_cov:
                continue

            hits.append(
                MMseqsHit(
                    query_id=query_id,
                    target_id=target_id,
                    target_header=target_raw,
                    identity=identity,
                    aln_len=aln_len,
                    query_cov=query_cov,
                    target_cov=target_cov,
                    qstart=int(parts[6]),
                    qend=int(parts[7]),
                    tstart=int(parts[8]),
                    tend=int(parts[9]),
                    evalue=float(parts[10]),
                    bits=float(parts[11]),
                )
            )

    return hits


def best_mmseqs_hits_by_family(
    hits: List[MMseqsHit],
) -> dict[str, MMseqsHit]:
    best = {}

    for hit in hits:
        current = best.get(hit.query_id)

        if current is None:
            best[hit.query_id] = hit
            continue

        better = False

        if hit.bits > current.bits:
            better = True

        elif hit.bits == current.bits:
            if hit.identity > current.identity:
                better = True

            elif (
                hit.identity == current.identity
                and hit.query_cov > current.query_cov
            ):
                better = True

            elif (
                hit.identity == current.identity
                and hit.query_cov == current.query_cov
                and hit.aln_len > current.aln_len
            ):
                better = True

        if better:
            best[hit.query_id] = hit

    return best
def reciprocal_best_hits(
    hits: List[MMseqsHit],
) -> List[MMseqsHit]:

    best = best_mmseqs_hits_by_family(hits)

    reciprocal = []

    for hit in hits:

        target_best = best.get(hit.target_id)

        if target_best is None:
            continue

        if (
            target_best.target_id == hit.query_id
        ):
            reciprocal.append(hit)

    return reciprocal

def infer_rescue_evidence(
    hits: List[MMseqsHit],
) -> dict[str, RescueEvidence]:

    reciprocal_hits = reciprocal_best_hits(
        hits
    )

    classified_hits = []

    for hit in reciprocal_hits:
        hom = homology_from_repeatmasker_header(
            hit.target_header
        )

        if hom is None:
            continue

        if hom.homology_order == "Unknown":
            continue

        classified_hits.append(hit)

    best_hits = best_mmseqs_hits_by_family(
        reciprocal_hits
    )

    rescue = {}

    for family_id, hit in best_hits.items():
        hom = homology_from_repeatmasker_header(
            hit.target_header
        )

        if hom is None:
            continue

        if hom.homology_order == "Unknown":
            continue

        if (
            hit.identity >= 0.80
            and hit.query_cov >= 0.80
            and hit.target_cov >= 0.80
        ):
            rescue_superfamily = hom.homology_superfamily
        else:
            rescue_superfamily = "Unknown"

        rescue[family_id] = RescueEvidence(
            family_id=family_id,
            rescue_class=hom.homology_class,
            rescue_order=hom.homology_order,
            rescue_superfamily=rescue_superfamily,
            rescue_identity=hit.identity,
            rescue_aln_len=hit.aln_len,
            rescue_query_cov=hit.query_cov,
            rescue_target_cov=hit.target_cov,
            rescue_bits=hit.bits,
            rescue_target=hit.target_header,
        )

    return rescue
