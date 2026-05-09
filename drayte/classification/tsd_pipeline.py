from __future__ import annotations

from pathlib import Path

from .ids import clean_family_id
from .tsd import TSDHit
from .tsd_from_alignment import detect_family_tsd_from_alignment


def detect_tsds_from_directory(
    align_dir: str | Path,
    suffix: str = ".fa",
    flank_size: int = 15,
    min_len: int = 2,
    max_len: int = 15,
    min_identity: float = 1.0,
    min_support_fraction: float = 0.5,
) -> list[TSDHit]:

    align_dir = Path(align_dir)

    if not align_dir.exists():
        return []

    hits: list[TSDHit] = []

    for fasta in sorted(align_dir.glob(f"*{suffix}")):
        family_id = clean_family_id(fasta.stem)

        hit = detect_family_tsd_from_alignment(
            family_id=family_id,
            alignment_fasta=fasta,
            flank_size=flank_size,
            min_len=min_len,
            max_len=max_len,
            min_identity=min_identity,
            min_support_fraction=min_support_fraction,
        )

        hits.append(hit)

    return hits
