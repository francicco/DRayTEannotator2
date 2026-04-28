from __future__ import annotations

from pathlib import Path

from drayte.utils.paths import stage_dir, ensure_dir
from drayte.structure import (
    run_ltr_module,
    run_tir_module,
    run_helitron_module,
    write_structure_library,
)
from drayte.structure.utils import write_candidates_tsv


def run(config, discovery_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "structure")

    logger.info("=" * 80)
    logger.info("STAGE: structure")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    genome_fa = config.outdir_path / "discovery" / "assemblies_dir" / f"{config.species}.fa"

    ltr_dir = ensure_dir(outdir / "ltr")
    tir_dir = ensure_dir(outdir / "tir")
    helitron_dir = ensure_dir(outdir / "helitron")

    all_candidates = []

    if config.extra.get("enable_ltr_module", True):
        all_candidates.extend(
            run_ltr_module(
                genome_fasta=genome_fa,
                outdir=ltr_dir,
                species=config.species,
                gt_bin=config.extra.get("genometools_bin", "gt"),
                logger=logger,
            )
        )

    if config.extra.get("enable_tir_module", False):
        all_candidates.extend(
            run_tir_module(
                genome_fasta=genome_fa,
                outdir=tir_dir,
                species=config.species,
                logger=logger,
            )
        )

    if config.extra.get("enable_helitron_module", False):
        all_candidates.extend(
            run_helitron_module(
                genome_fasta=genome_fa,
                outdir=helitron_dir,
                species=config.species,
                logger=logger,
            )
        )

    manifest_tsv = outdir / "structure_candidates.tsv"
    write_candidates_tsv(all_candidates, manifest_tsv)

    structure_library = outdir / f"{config.species}_structure_candidates.fa"
    write_structure_library(all_candidates, structure_library, logger)

    result = {
        "stage": "structure",
        "outdir": str(outdir),
        "structure_library": str(structure_library),
        "structure_manifest": str(manifest_tsv),
        "n_structure_candidates": len(all_candidates),
    }

    logger.info("Structure stage completed")
    return result
