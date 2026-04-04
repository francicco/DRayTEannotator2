from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class StructureCandidate:
    candidate_id: str
    module: str
    contig: str
    start: int
    end: int
    strand: str
    length: int
    fasta_path: Path
