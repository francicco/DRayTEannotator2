from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from drayte.utils.paths import stage_dir, ensure_dir


def final_annotation_exists(outdir: Path, taxon: str) -> bool:
    required = [
        outdir / f"{taxon}.fa.out",
        outdir / f"{taxon}.fa.tbl",
        outdir / f"{taxon}.fa.masked",
        outdir / f"{taxon}.fa.cat.gz",
        outdir / f"{taxon}.fa.align",
    ]

    return all(f.exists() and f.stat().st_size > 0 for f in required)


def move_repeatmasker_outputs(genome_fa: Path, outdir: Path, logger) -> None:
    parent = genome_fa.parent
    genome_name = genome_fa.name

    for suffix in [".out", ".tbl", ".gff", ".masked", ".cat.gz", ".align"]:
        src = parent / f"{genome_name}{suffix}"
        dst = outdir / f"{genome_name}{suffix}"

        if src.exists():
            if dst.exists():
                dst.unlink()
            logger.info("Moving %s -> %s", src, dst)
            shutil.move(str(src), str(dst))


def run(config, curation_result: dict, logger) -> dict:
    outdir = ensure_dir(stage_dir(config.outdir_path, "final_annotation"))

    logger.info("=" * 80)
    logger.info("STAGE: final_annotation")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    taxon = config.species
    genome_fa = config.outdir_path / "discovery" / "assemblies_dir" / f"{taxon}.fa"
    final_library = Path(curation_result["final_library"])

    tbl = outdir / f"{taxon}.fa.tbl"
    out = outdir / f"{taxon}.fa.out"
    gff = outdir / f"{taxon}.fa.gff"
    masked = outdir / f"{taxon}.fa.masked"
    cat_gz = outdir / f"{taxon}.fa.cat.gz"
    align = outdir / f"{taxon}.fa.align"

    if final_annotation_exists(outdir, taxon):
        logger.info("Final RepeatMasker outputs already exist; skipping")
    else:
        cmd = [
            str(config.extra.get("repeatmasker_bin", "RepeatMasker")),
            "-s",
            "-xsmall",
            "-gc",
            "-gff",
            "-norna",
            "-a",
            "-pa", str(config.threads),
            "-lib", str(final_library),
            str(genome_fa),
        ]

        logger.info("Running final RepeatMasker: %s", " ".join(cmd))
        subprocess.run(cmd, check=True, cwd=str(genome_fa.parent))
        move_repeatmasker_outputs(genome_fa, outdir, logger)

    return {
        "stage": "final_annotation",
        "outdir": str(outdir),
        "repeatmasker_tbl": str(tbl),
        "repeatmasker_out": str(out),
        "repeatmasker_gff": str(gff),
        "masked_genome": str(masked),
        "cat_gz": str(cat_gz),
        "align": str(align),
    }
