from pathlib import Path

from drayte.utils.paths import stage_dir
from drayte.utils.subprocess import run_command


def run(config, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "discovery")

    logger.info("Starting discovery stage")

    cmd = [
        "python",
        config.tools.step1_repmodannotation,
        "--genome",
        str(config.genome_path.resolve()),
        "--species",
        config.species,
        "--threads",
        str(config.threads),
        "--outdir",
        str(outdir),
    ]

    run_command(cmd, logger)

    result = {
        "stage": "discovery",
        "outdir": str(outdir),
        "genome": str(config.genome_path.resolve()),
        "raw_library": str(outdir / "consensi.fa.classified"),
    }

    logger.info("Discovery stage completed")
    return result
