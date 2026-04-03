from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Dict

from Bio import SeqIO

from .metadata import FamilyMeta


def run_teaid(
    teaid_bin: str,
    query_fa: Path,
    genome_fa: Path,
    outdir: Path,
    logger,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        teaid_bin,
        "-q", str(query_fa),
        "-g", str(genome_fa),
        "-T",
        "-o", str(outdir),
    ]
    logger.info("Running TE-Aid on %s", query_fa.name)
    subprocess.run(cmd, check=True)

    pdf = Path(f"{query_fa}.c2g.pdf")
    if pdf.exists():
        pdf.rename(outdir / f"{query_fa.stem}.c2g.pdf")


def count_full_length_hits(blastn_out: Path, consensus_len: int) -> int:
    if not blastn_out.exists():
        return 0

    min_len = consensus_len * 0.9
    count = 0
    with open(blastn_out) as handle:
        next(handle, None)
        for line in handle:
            row = line.strip().split()
            if len(row) < 8:
                continue
            s = int(row[6])
            e = int(row[7])
            alen = abs(e - s) + 1
            if alen > min_len:
                count += 1
    return count


def append_teaid_rows(
    metadata: Dict[str, FamilyMeta],
    category_of: Dict[str, str],
    aidout: Path,
    outfile: Path,
) -> Path:
    with open(outfile, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "RM_ID", "Short_ID", "Class", "Family", "Modified_ID",
            "Consensus_length", "90percent_consensus",
            "N_ORFS", "ORF1_type", "ORF1_length",
            "ORF2_type", "ORF2_length", "ORF3_type", "ORF3_length",
        ])

        for short_id, meta in sorted(metadata.items()):
            category = category_of[short_id]
            rep_mod = aidout / category / f"{short_id}_rep_mod.fa"
            modified_id = short_id
            if rep_mod.exists():
                with open(rep_mod) as f:
                    header = f.readline().strip()
                modified_id = header[1:] if header.startswith(">") else short_id

            blastn_out = aidout / category / f"{short_id}_rep.fa.genome.blastn.out"
            fullcount = count_full_length_hits(blastn_out, meta.consensus_len)

            writer.writerow([
                meta.original_header, short_id, meta.te_class, meta.te_family, modified_id,
                meta.consensus_len, fullcount
