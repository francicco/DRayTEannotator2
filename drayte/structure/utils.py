from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta

from .models import StructureCandidate


def write_candidates_tsv(candidates: list[StructureCandidate], outfile: Path) -> None:
    with open(outfile, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["candidate_id", "module", "contig", "start", "end", "strand", "length", "fasta_path"])
        for c in candidates:
            writer.writerow([
                c.candidate_id,
                c.module,
                c.contig,
                c.start,
                c.end,
                c.strand,
                c.length,
                str(c.fasta_path),
            ])


def extract_region_fasta(
    genome_fasta: Path,
    contig: str,
    start_1based: int,
    end_1based: int,
    strand: str,
    out_fasta: Path,
    seq_id: str,
) -> Path:
    genome = Fasta(str(genome_fasta), as_raw=True, sequence_always_upper=True)

    start0 = max(0, start_1based - 1)
    end0 = end_1based
    seq = str(genome[contig][start0:end0])

    if strand == "-":
        seq = str(Seq(seq).reverse_complement())

    rec = SeqRecord(Seq(seq), id=seq_id, description="")
    with open(out_fasta, "w") as handle:
        SeqIO.write([rec], handle, "fasta")

    return out_fasta


def concatenate_fastas(fastas: Iterable[Path], outfile: Path) -> int:
    count = 0
    with open(outfile, "w") as out_handle:
        for fasta in fastas:
            if not fasta.exists() or fasta.stat().st_size == 0:
                continue
            for rec in SeqIO.parse(str(fasta), "fasta"):
                SeqIO.write(rec, out_handle, "fasta")
                count += 1
    return count
