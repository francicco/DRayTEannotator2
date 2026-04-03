from drayte.utils.paths import stage_dir
from drayte.utils.subprocess import run_command


def run(config, reclassify_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "curation")

    logger.info("Starting curation stage")

    cmd = [
        "bash",
        config.tools.step4_tecurate,
        reclassify_result["classified_library"],
        str(config.genome_path.resolve()),
        str(outdir),
    ]

    run_command(cmd, logger)

    result = {
        "stage": "curation",
        "outdir": str(outdir),
        "curated_library": str(outdir / "Final.RepeatModeler.Lib.fa"),
    }

    logger.info("Curation stage completed")
    return result
