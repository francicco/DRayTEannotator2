from pathlib import Path

from drayte.utils.paths import stage_dir
from drayte.step1_repmodannotation import run_step1


def run(config, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "discovery")

    logger.info("=" * 80)
    logger.info("STAGE: discovery")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    step1_result = run_step1(
        genome=config.genome_path.resolve(),
        outdir=outdir,
        species=config.species,
        threads=config.threads,
        repeatmodeler_dir=Path(config.extra["repeatmodeler_dir"]),
        repeatscout_dir=Path(config.extra["repeatscout_dir"]),
        repeatmasker_bin=config.extra.get("repeatmasker_bin", "RepeatMasker"),
        logger=logger,
    )

    result = {
        "stage": "discovery",
        "outdir": str(outdir),
        "genome": str(config.genome_path.resolve()),
        "raw_library": step1_result["raw_library"],
        "normalized_library": step1_result["normalized_library"],
        "repeatmasker_dir": step1_result["rmasker_dir"],
    }

    logger.info("Discovery stage completed")
    return result
