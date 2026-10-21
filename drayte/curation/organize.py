from __future__ import annotations

import csv
import shutil
from pathlib import Path

from drayte.utils.names import safe_filename

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

    extension_dirs: dict[str, list[Path]] = {}
    if extensionwork_dir.exists():
        for candidate in extensionwork_dir.iterdir():
            if not candidate.is_dir():
                continue
            canonical_name = candidate.name.split("#", 1)[0]
            extension_dirs.setdefault(canonical_name, []).append(candidate)

    matched = 0
    missing = 0
    ambiguous = 0

    with open(family_table_tsv) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            consname = row["name"]
            cons_short = consname
            group = classify_group(row)

            canonical_name = safe_filename(cons_short).split("#", 1)[0]
            candidates = extension_dirs.get(canonical_name, [])

            if not candidates:
                missing += 1
                logger.warning(
                    "Missing extension directory for %s under %s",
                    consname,
                    extensionwork_dir,
                )
                continue

            if len(candidates) > 1:
                ambiguous += 1
                logger.error(
                    "Ambiguous extension directories for %s: %s",
                    consname,
                    ", ".join(str(path) for path in sorted(candidates)),
                )
                continue

            src_dir = candidates[0]
            source_prefix = src_dir.name
            matched += 1

            safe_consname = safe_filename(consname)

            mapping = {
                f"{source_prefix}_rep.fa": te_dirs[group] / f"{safe_consname}_rep.fa",
                f"{source_prefix}_MSA_extended.fa": te_dirs[group] / f"{safe_consname}_MSA_extended.fa",
                f"{source_prefix}.png": te_dirs[group] / f"{safe_consname}.png",
            }

            found_any = False
            
            for src_name, dst in mapping.items():
                src = src_dir / src_name
                if src.exists():
                    found_any = True
                    if not dst.exists():
                        shutil.copy2(src, dst)
            
            if not found_any:
                logger.debug(
                    "No extension artifacts for %s (likely low-copy or skipped)",
                    consname,
                )

    logger.info(
        "Extension artifact directories: %d matched, %d missing, %d ambiguous",
        matched,
        missing,
        ambiguous,
    )

    if ambiguous:
        raise RuntimeError(
            f"Found {ambiguous} ambiguous extension directory matches; "
            "refusing to choose an arbitrary directory"
        )

    return te_dirs
