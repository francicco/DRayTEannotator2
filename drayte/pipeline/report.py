from drayte.utils.paths import stage_dir


def run(config, curation_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "report")

    logger.info("Starting report stage")
    logger.info("Report stage is currently a placeholder")

    result = {
        "stage": "report",
        "outdir": str(outdir),
        "curated_library": curation_result["curated_library"],
    }

    logger.info("Report stage completed")
    return result
