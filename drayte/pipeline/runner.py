import argparse
import json
from pathlib import Path

from drayte.pipeline.config import load_config
from drayte.pipeline import discovery, extension, reclassify, curation, final_annotation, report
from drayte.structure import heliano
from drayte.utils.logging import setup_logger
from drayte.utils.paths import ensure_dir
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

    final_annotation_result = final_annotation.run(config, curation_result, logger)
    write_manifest(config.outdir_path / "final_annotation", "final_annotation", final_annotation_result)

    summary_result = run_summary_files(
        config=config,
        final_annotation_result=final_annotation_result,
        logger=logger,
    )
    write_manifest(config.outdir_path / "summaryFiles", "summaryFiles", summary_result)

    report_result = report.run(config, curation_result, logger)
    write_manifest(config.outdir_path / "report", "report", report_result)

    logger.info("Pipeline completed successfully")

if __name__ == "__main__":
    main()
