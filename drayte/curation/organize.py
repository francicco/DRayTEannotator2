from __future__ import annotations

import csv
import shutil
from pathlib import Path

from drayte.extension.extract_align import safe_filename

TELIST = ["LINE", "SINE", "LTR", "RC", "DNA", "NOHIT"]


def ensure_te_dirs(base_dir: Path) -> dict[str, Path]:
    base_dir.mkdir(parents=True, exist_ok=True)
    dirs = {}
    for te in TELIST:
        d = base_dir / te
        d.mkdir(parents=True, exist_ok=True)
        dirs[te] = d
    return dirs


def classify_group(row: dict) -> str:
    top_orf_class = row["top_orf_class"]
    if top_orf_class in {"LINE", "SINE", "LTR", "RC", "DNA"}:
        return top_orf_class
    return "NOHIT"


def copy_extension_artifacts(
    family_table_tsv: Path,
    extensionwork_dir: Path,
    outdir: Path,
    logger,
) -> dict[str, Path]:
    te_dirs = ensure_te_dirs(outdir)

    with open(family_table_tsv) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            consname = row["name"]
            cons_short = consname
            group = classify_group(row)

            src_dir = extensionwork_dir / safe_filename(cons_short)
            if not src_dir.exists():
                logger.warning("Missing extension directory for %s: %s", consname, src_dir)
                continue

            mapping = {
                f"{cons_short}_rep.fa": te_dirs[group] / f"{consname}_rep.fa",
                f"{cons_short}_MSA_extended.fa": te_dirs[group] / f"{consname}_MSA_extended.fa",
                f"{cons_short}.png": te_dirs[group] / f"{consname}.png",
            }

            for src_name, dst in mapping.items():
                src = src_dir / src_name
                if src.exists():
                    if not dst.exists():
                        shutil.copy2(src, dst)
                else:
                    logger.warning("Missing expected extension artifact: %s", src)

    return te_dirs
