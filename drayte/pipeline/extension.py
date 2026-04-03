from drayte.utils.paths import stage_dir
from drayte.utils.subprocess import run_command


def run(config, discovery_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "extension")

    logger.info("Starting extension stage")

    cmd = [
        "bash",
        config.tools.step2_extendedalign,
        discovery_result["raw_library"],
        str(config.genome_path.resolve()),
        str(outdir),
    ]

    run_command(cmd, logger)

    result = {
        "stage": "extension",
        "outdir": str(outdir),
        "extended_library": str(outdir / "extended_consensi.fa"),
    }

    logger.info("Extension stage completed")
    return result
