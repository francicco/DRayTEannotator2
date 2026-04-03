from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO


def build_family_table(
    library_fasta: Path,
    orf_result: dict,
    homology_result: dict,
    outdir: Path,
    logger,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    table = outdir / "families.tsv"

    records = list(SeqIO.parse(str(library_fasta), "fasta"))

    with open(table, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "family_id",
            "consensus_len",
            "orf_fasta",
            "diamond_tsv",
        ])
        for rec in records:
            writer.writerow([
                rec.id,
                len(rec.seq),
                orf_result.get("orf_fasta"),
                homology_result.get("diamond_tsv"),
            ])

    logger.info("Wrote family table: %s", table)
    return table
