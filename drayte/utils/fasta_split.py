from __future__ import annotations

from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def split_fasta_balanced_by_length(
    input_fasta: str | Path,
    outdir: str | Path,
    n_chunks: int,
    prefix: str = "chunk",
) -> list[Path]:
    """
    Split a FASTA file into approximately balanced chunks by total sequence length.

    Sequences are sorted by length descending and greedily assigned to the
    currently smallest chunk. This is simple bin packing and generally gives
    better wall-time balance than splitting by number of records.
    """
    input_fasta = Path(input_fasta)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if n_chunks < 1:
        raise ValueError("n_chunks must be >= 1")

    records = list(SeqIO.parse(str(input_fasta), "fasta"))

    if not records:
        return []

    records.sort(
        key=lambda rec: len(rec.seq),
        reverse=True,
    )

    n_chunks = min(n_chunks, len(records))

    bins: list[list[SeqRecord]] = [
        [] for _ in range(n_chunks)
    ]

    bin_sizes = [0] * n_chunks

    for rec in records:
        idx = bin_sizes.index(min(bin_sizes))
        bins[idx].append(rec)
        bin_sizes[idx] += len(rec.seq)

    chunk_paths = []

    width = len(str(n_chunks))

    for i, chunk_records in enumerate(bins, start=1):
        chunk_path = outdir / f"{prefix}_{i:0{width}d}.fa"

        with open(chunk_path, "w") as out:
            SeqIO.write(chunk_records, out, "fasta")

        chunk_paths.append(chunk_path)

    report_path = outdir / f"{prefix}.sizes.tsv"

    with open(report_path, "w") as out:
        out.write("chunk\trecords\tbp\n")

        for path, chunk_records, size in zip(
            chunk_paths,
            bins,
            bin_sizes,
        ):
            out.write(
                f"{path.name}\t{len(chunk_records)}\t{size}\n"
            )

    return chunk_paths
