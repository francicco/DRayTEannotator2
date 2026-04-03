from pathlib import Path
from typing import Iterator
from Bio import SeqIO


def filter_fasta_by_length(
    input_fasta: str | Path,
    output_fasta: str | Path,
    min_length: int = 100,
) -> int:
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)

    kept = []
    for record in SeqIO.parse(str(input_fasta), "fasta"):
        if len(record.seq) >= min_length:
            kept.append(record)

    with open(output_fasta, "w") as handle:
        SeqIO.write(kept, handle, "fasta")

    return len(kept)


def rename_repeatmodeler_headers(
    input_fasta: str | Path,
    output_fasta: str | Path,
    species: str,
) -> int:
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)

    records = []
    for record in SeqIO.parse(str(input_fasta), "fasta"):
        header = record.id.split()[0]
        header = header.replace("#", "__")
        header = header.replace("/", "___")
        header = header.replace("rnd-", f"{species}-rnd-")
        record.id = header
        record.description = header
        records.append(record)

    with open(output_fasta, "w") as handle:
        SeqIO.write(records, handle, "fasta")

    return len(records)
