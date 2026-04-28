from __future__ import annotations

from pathlib import Path

from .models import StructureCandidate
from .utils import concatenate_fastas


def write_structure_library(
    candidates: list[StructureCandidate],
    outfile: Path,
    logger,
) -> int:
    count = concatenate_fastas([c.fasta_path for c in candidates], outfile)
    logger.info("Wrote structure library with %d candidates: %s", count, outfile)
    return count
