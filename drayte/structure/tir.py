from __future__ import annotations

from pathlib import Path

from .models import StructureCandidate


def run_tir_module(
    genome_fasta: Path,
    outdir: Path,
    species: str,
    logger,
) -> list[StructureCandidate]:
    outdir.mkdir(parents=True, exist_ok=True)
    logger.info("TIR module is not implemented yet; returning zero candidates")
    return []
