from drayte.utils.paths import stage_dir
from drayte.utils.subprocess import run_command


def run(config, extension_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "reclassify")

    logger.info("Starting reclassify stage")

    cmd = [
        "bash",
        config.tools.step3_repeatclassifier,
        extension_result["extended_library"],
        str(outdir),
    ]

    run_command(cmd, logger)

    result = {
        "stage": "reclassify",
        "outdir": str(outdir),
        "classified_library": str(outdir / "extended_consensi.fa.classified"),
    }

    logger.info("Reclassify stage completed")
    return result
