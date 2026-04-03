from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict

from .metadata import FamilyMeta


CLASS_MAP = {
    "LINE": "LINE",
    "SINE": "SINE",
    "LTR": "LTR",
    "RC": "RC",
    "DNA": "DNA",
    "Unknown": "NOHIT",
}


def short_id_to_mod_id(short_id: str, te_class: str, te_family: str) -> str:
    consnamemod = short_id.replace("-rnd-", ".").replace("_family-", ".")
    return f"{consnamemod}#{te_class}/{te_family}"


def organize_family_outputs(
    metadata: Dict[str, FamilyMeta],
    extension_workdir: str | Path,
    aidout: str | Path,
    logger,
) -> Dict[str, str]:
    extension_workdir = Path(extension_workdir)
    aidout = Path(aidout)
    for cat in ["LINE", "SINE", "LTR", "RC", "DNA", "NOHIT"]:
        (aidout / cat).mkdir(parents=True, exist_ok=True)

    category_of: Dict[str, str] = {}

    for short_id, meta in metadata.items():
        category = CLASS_MAP.get(meta.te_class, "NOHIT")
        category_of[short_id] = category

        src_dir = extension_workdir / short_id
        dst_dir = aidout / category

        rep = src_dir / f"{short_id}_rep.fa"
        msa = src_dir / f"{short_id}_MSA_extended.fa"
        png = src_dir / f"{short_id}.png"

        if rep.exists():
            shutil.copy2(rep, dst_dir / f"{short_id}_rep.fa")
        if msa.exists():
            shutil.copy2(msa, dst_dir / f"{short_id}_MSA_extended.fa")
        if png.exists():
            shutil.copy2(png, dst_dir / f"{short_id}.png")

        rep_mod = dst_dir / f"{short_id}_rep_mod.fa"
        if (dst_dir / f"{short_id}_rep.fa").exists() and not rep_mod.exists():
            new_header = short_id_to_mod_id(short_id, meta.te_class, meta.te_family)
            with open(dst_dir / f"{short_id}_rep.fa") as fin, open(rep_mod, "w") as fout:
                first = True
                for line in fin:
                    if first and line.startswith(">"):
                        fout.write(f">{new_header}\n")
                        first = False
                    else:
                        fout.write(line)

    logger.info("Organized TE-Aid input files by class")
    return category_of
