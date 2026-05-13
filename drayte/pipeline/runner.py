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
from drayte.refinement.refine_repeatmasker import refine_repeatmasker
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

	# Run RepeatModeler
    discovery_result = discovery.run(config, logger)
    write_manifest(config.outdir_path / "discovery", "discovery", discovery_result)

	# Run Extention - Part 1
    extension_result = extension.run(config, discovery_result, logger)
    write_manifest(config.outdir_path / "extension", "extension", extension_result)

	# Run Re-classify RepeatModeler libray using extentions - Part 2
    reclassify_result = reclassify.run(config, extension_result, logger)
    write_manifest(config.outdir_path / "reclassify", "reclassify", reclassify_result)

	# Run curation - Part 3
    curation_result = curation.run(config, reclassify_result, logger)
    write_manifest(config.outdir_path / "curation", "curation", curation_result)

	# Optional - Run HELIANO to identify full length Helitron 1 and 2
    heliano_result = heliano.run(config, curation_result, logger)
    write_manifest(config.outdir_path / "heliano", "heliano", heliano_result)

    if heliano_result.get("heliano_unique_library"):
        curation_result["heliano_library"] = heliano_result["heliano_unique_library"]

	# Run non-autonomous TEs - New Implementation
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

	# Rerun RepeatMasker with the liberary just classified further
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

    refinement_outdir = config.outdir_path / "Defragmenting_RModAnnotation"

    refined_tsv = refinement_outdir / f"{config.species}.annotation_refinement.tsv"
    filtered_tsv = refinement_outdir / f"{config.species}.filteredRepeats.tsv"
    stats_tsv = refinement_outdir / f"{config.species}.annotation_refinement.stats.tsv"

    if (
        refined_tsv.exists() and refined_tsv.stat().st_size > 0
        and filtered_tsv.exists() and filtered_tsv.stat().st_size > 0
        and stats_tsv.exists() and stats_tsv.stat().st_size > 0
    ):
        logger.info("=" * 80)
        logger.info("STAGE: annotation_refinement")
        logger.info("Output directory: %s", refinement_outdir)
        logger.info("=" * 80)
        logger.info("TE-refine outputs already exist; skipping")

        refinement_result = {
            "stage": "annotation_refinement",
            "outdir": str(refinement_outdir),
            "refined_tsv": str(refined_tsv),
            "filtered_tsv": str(filtered_tsv),
            "stats_tsv": str(stats_tsv),
            "refined_gff": str(refinement_outdir / f"{config.species}.refinedRepeats.gff3"),
            "refined_bed": str(refinement_outdir / f"{config.species}.refinedRepeats.bed"),
            "filtered_gff": str(refinement_outdir / f"{config.species}.filteredRepeats.gff3"),
            "filtered_bed": str(refinement_outdir / f"{config.species}.filteredRepeats.bed"),
            "nested_gff": str(refinement_outdir / f"{config.species}.nestedRepeats.gff3"),
        }
    else:
        refinement_result = refine_repeatmasker(
            rmout=Path(final_annotation_result["repeatmasker_out"]),
            species=config.species,
            outdir=refinement_outdir,
            max_gap=int(config.extra.get("refinement_max_gap", 150)),
            include_nested=bool(config.extra.get("refinement_include_nested", False)),
            mode=config.extra.get("refinement_mode", "loose"),
            max_consensus_overlap=int(config.extra.get("refinement_max_consensus_overlap", 50)),
            max_consensus_jump=int(config.extra.get("refinement_max_consensus_jump", 10000)),
            resolve_overlaps=bool(config.extra.get("refinement_resolve_overlaps", True)),
            overlap_score=config.extra.get("refinement_overlap_score", "longest_lowdiv"),
            logger=logger,
        )

    write_manifest(
        refinement_outdir,
        "annotation_refinement",
        refinement_result,
    )

    refinement_outdir = config.outdir_path / "annotation_refinement"

    summary_no_filtering_outdir = (
        config.outdir_path.parent / f"{config.species}.annotation_refinement_noFiltering"
    )

    summary_no_filter_result = run_summary_files(
        refined_tsv=config.outdir_path / "annotation_refinement" / f"{config.species}.annotation_refinement.tsv",
        genome=Path(config.genome),
        species=config.species,
        outdir=config.outdir_path / f"{config.species}.annotation_refinement_noFiltered",
        max_merge_gap=int(config.extra.get("summary_max_merge_gap", 100)),
        min_nested_overlap_fraction=float(config.extra.get("summary_min_nested_overlap_fraction", 0.80)),
        logger=logger,
    )
    
    summary_filtered_result = run_summary_files(
        refined_tsv=config.outdir_path / "annotation_refinement" / f"{config.species}.filteredRepeats.tsv",
        genome=Path(config.genome),
        species=config.species,
        outdir=config.outdir_path / f"{config.species}.annotation_refinement_Filtered",
        max_merge_gap=int(config.extra.get("summary_max_merge_gap", 100)),
        min_nested_overlap_fraction=float(config.extra.get("summary_min_nested_overlap_fraction", 0.80)),
        logger=logger,
    )

    write_manifest(
        summary_no_filtering_outdir,
        "TE-summary.noFiltering",
        summary_no_filter_result,
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
