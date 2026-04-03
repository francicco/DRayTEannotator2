from __future__ import annotations

import shutil
from pathlib import Path

from Bio import SeqIO


def strip_classification_from_headers(input_fasta: Path, output_fasta: Path, logger) -> Path:
    records = []
    for record in SeqIO.parse(str(input_fasta), "fasta"):
        record.id = record.id.split("#")[0]
        record.description = record.id
        records.append(record)

    with open(output_fasta, "w") as handle:
        SeqIO.write(records, handle, "fasta")

    logger.info("Wrote cleaned library FASTA: %s", output_fasta)
    return output_fasta


def prepare_curation_inputs(classified_library: Path, outdir: Path, species: str, logger) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)

    copied = outdir / f"{species}_extended_rep.fa.classified"
    if not copied.exists():
        shutil.copy2(classified_library, copied)

    clean_library = outdir / f"{species}_extended_rep.clean.fa"
    if not clean_library.exists():
        strip_classification_from_headers(copied, clean_library, logger)

    return {
        "classified_library": copied,
        "clean_library": clean_library,
    }
