import argparse
import json
from pathlib import Path

from drayte.pipeline.config import load_config
from drayte.structure import heliano
from drayte.utils.logging import setup_logger
from drayte.utils.paths import ensure_dir
from drayte.reporting import family_inspection
from drayte.pipeline import (
    discovery,
    extension,
    reclassify,
    curation,
    classification,
    final_annotation,
    report,
    annotation_refinement,
)
from drayte.reporting.SummaryFilesGen import run_summary_files

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="DRayTE pipeline runner")
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    return parser.parse_args()


def write_manifest(outdir: Path, name: str, data: dict) -> None:
    outpath = outdir / f"{name}.manifest.json"
    with open(outpath, "w") as handle:
        json.dump(data, handle, indent=2)


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    ensure_dir(config.outdir_path)
    logger = setup_logger(log_file=str(config.outdir_path / "drayte.log"))

    logger.info("Loaded config from %s", args.config)
    logger.info("Running pipeline for species=%s", config.species)

    discovery_result = discovery.run(config, logger)
    write_manifest(config.outdir_path / "discovery", "discovery", discovery_result)

    extension_result = extension.run(config, discovery_result, logger)
    write_manifest(config.outdir_path / "extension", "extension", extension_result)

    reclassify_result = reclassify.run(config, extension_result, logger)
    write_manifest(config.outdir_path / "reclassify", "reclassify", reclassify_result)

    curation_result = curation.run(config, reclassify_result, logger)
    write_manifest(config.outdir_path / "curation", "curation", curation_result)

    heliano_result = heliano.run(config, curation_result, logger)
    write_manifest(config.outdir_path / "heliano", "heliano", heliano_result)

    if heliano_result.get("heliano_unique_library"):
        curation_result["heliano_library"] = heliano_result["heliano_unique_library"]

    classification_result = classification.run(
        config,
        curation_result,
        logger,
        stage_name="classification",
        final_mode=True,
    )
    write_manifest(
        config.outdir_path / "classification",
        "classification",
        classification_result,
    )

    final_annotation_input = dict(curation_result)
    final_annotation_input["classified_library"] = classification_result["classified_library"]

    final_annotation_result = final_annotation.run(
        config,
        final_annotation_input,
        logger,
    )
    write_manifest(
        config.outdir_path / "final_annotation",
        "final_annotation",
        final_annotation_result,
    )

    refinement_result = annotation_refinement.run(
        config,
        final_annotation_result,
        logger,
    )
    write_manifest(
        config.outdir_path / "annotation_refinement",
        "annotation_refinement",
        refinement_result,
    )

    family_inspection_result = family_inspection.run(
        config,
        refinement_result,
        logger,
    )
    write_manifest(
        config.outdir_path / "family_inspection",
        "family_inspection",
        family_inspection_result,
    )

    logger.info("Pipeline completed successfully")

if __name__ == "__main__":
    main()
