from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO


@dataclass
class FamilyMeta:
    original_header: str
    rm_name: str
    te_class: str
    te_family: str
    consensus_len: int = 0
    n_orfs: int = 0
    orf1_type: str = "NOHIT"
    orf1_len: int = 0
    orf2_type: str = "NOHIT"
    orf2_len: int = 0
    orf3_type: str = "NOHIT"
    orf3_len: int = 0


def _normalise_header(header: str) -> tuple[str, str, str]:
    clean = header.replace("#", "-#")
    if "#Unknown" in clean:
        clean = clean.replace("#Unknown", "#Unknown/Unknown")
    if "#Satellite" in clean:
        clean = clean.replace("#Satellite", "#Satellite/Satellite")
    if "#LTR " in clean:
        clean = clean.replace("#LTR ", "#LTR/Unknown")
    if "#DNA " in clean:
        clean = clean.replace("#DNA ", "#DNA/Unknown")
    if "#tRNA " in clean:
        clean = clean.replace("#tRNA ", "#tRNA/Nothing")
    if "#LINE " in clean:
        clean = clean.replace("#LINE ", "#LINE/Unknown")

    rm_name = clean.split("#")[0]
    class_family = clean.split("#", 1)[1] if "#" in clean else "Unknown/Unknown"
    if "/" in class_family:
        te_class, te_family = class_family.split("/", 1)
    else:
        te_class, te_family = class_family, "Unknown"

    return rm_name, te_class, te_family


def load_family_metadata(classified_fasta: str | Path) -> Dict[str, FamilyMeta]:
    classified_fasta = Path(classified_fasta)
    meta: Dict[str, FamilyMeta] = {}

    for record in SeqIO.parse(str(classified_fasta), "fasta"):
        original_header = record.id
        rm_name, te_class, te_family = _normalise_header(original_header)
        meta[rm_name] = FamilyMeta(
            original_header=original_header,
            rm_name=rm_name,
            te_class=te_class,
            te_family=te_family,
            consensus_len=len(record.seq),
        )
    return meta


def write_family_table(metadata: Dict[str, FamilyMeta], outfile: str | Path) -> Path:
    outfile = Path(outfile)
    with open(outfile, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "RM_ID", "Short_ID", "Class", "Family",
            "Consensus_length", "N_ORFS",
            "ORF1_type", "ORF1_length",
            "ORF2_type", "ORF2_length",
            "ORF3_type", "ORF3_length",
        ])
        for rm_name, m in sorted(metadata.items()):
            writer.writerow([
                m.original_header,
                rm_name,
                m.te_class,
                m.te_family,
                m.consensus_len,
                m.n_orfs,
                m.orf1_type,
                m.orf1_len,
                m.orf2_type,
                m.orf2_len,
                m.orf3_type,
                m.orf3_len,
            ])
    return outfile
