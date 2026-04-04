from __future__ import annotations

from pathlib import Path

from drayte.utils.paths import stage_dir, ensure_dir
from drayte.curation import (
    build_family_table_from_classified_library,
    copy_extension_artifacts,
    orient_group_files,
)


def run(config, reclassify_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "curation")

    logger.info("=" * 80)
    logger.info("STAGE: curation")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    classified_library = Path(reclassify_result["classified_library"])
    extension_dir = config.outdir_path / "extension"
    extensionwork_dir = extension_dir / "extensionwork"

    prioritize_dir = ensure_dir(outdir / "prioritize")
    te_aid_dir = ensure_dir(outdir / "te-aid")

    family_table = build_family_table_from_classified_library(
        classified_library=classified_library,
        outdir=prioritize_dir,
        diamond_bin=config.extra.get("diamond_bin", "diamond"),
        repeatpeps_db_dir=Path(config.extra["repeatpeps_db_dir"]),
        min_orf=int(config.extra.get("min_orf", 500)),
        threads=config.threads,
        logger=logger,
    )

    copy_extension_artifacts(
        family_table_tsv=family_table,
        extensionwork_dir=extensionwork_dir,
        outdir=te_aid_dir,
        logger=logger,
    )

    orient_group_files(
        family_table_tsv=family_table,
        te_aid_dir=te_aid_dir,
        diamond_bin=config.extra.get("diamond_bin", "diamond"),
        repeatpeps_db_dir=Path(config.extra["repeatpeps_db_dir"]),
        threads=config.threads,
        logger=logger,
    )

    result = {
        "stage": "curation",
        "outdir": str(outdir),
        "family_table": str(family_table),
        "te_aid_dir": str(te_aid_dir),
    }

    logger.info("Curation stage completed")
    return result
