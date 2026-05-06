from pathlib import Path
from typing import Dict, List

from Bio import SeqIO

from drayte.curation.orf_compat import (
    OrfCall,
    find_orfs_getorf_compatible,
    write_orfs_fasta,
)


def extract_orf_calls_by_family(
    input_fasta: str | Path,
    minsize_nt: int = 500,
    include_reverse: bool = True,
) -> Dict[str, List[OrfCall]]:
    calls_by_family: Dict[str, List[OrfCall]] = {}

    for rec in SeqIO.parse(str(input_fasta), "fasta"):
        calls_by_family[rec.id] = find_orfs_getorf_compatible(
            seq_record=rec,
            minsize_nt=minsize_nt,
            include_reverse=include_reverse,
        )

    return calls_by_family


def summarize_orfs(
    calls_by_family: Dict[str, List[OrfCall]],
) -> Dict[str, dict]:
    summary = {}

    for family_id, calls in calls_by_family.items():
        nt_lengths = [c.nt_length for c in calls]

        summary[family_id] = {
            "orf_count": len(calls),
            "orf_max_len": max(nt_lengths) if nt_lengths else 0,
        }

    return summary


def write_translated_orfs(
    input_fasta: str | Path,
    output_fasta: str | Path,
    minsize_nt: int = 500,
    include_reverse: bool = True,
) -> int:
    return write_orfs_fasta(
        input_fasta=Path(input_fasta),
        output_fasta=Path(output_fasta),
        minsize_nt=minsize_nt,
        include_reverse=include_reverse,
    )
