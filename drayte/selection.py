from __future__ import annotations

import shutil
from pathlib import Path


def select_representatives(
    family_table_tsv: Path,
    library_fasta: Path,
    outdir: Path,
    logger,
) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)

    curated_library = outdir / "Final.RepeatModeler.Lib.fa"
    selected_table = outdir / "selected_families.tsv"

    if not curated_library.exists():
        shutil.copy2(library_fasta, curated_library)

    if not selected_table.exists():
        shutil.copy2(family_table_tsv, selected_table)

    logger.info("Selected representative library: %s", curated_library)

    return {
        "curated_library": str(curated_library),
        "selected_table": str(selected_table),
    }
