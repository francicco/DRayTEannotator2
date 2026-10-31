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
from drayte.refinement.refine_repeatmasker import (
    REFINEMENT_ALGORITHM_VERSION,
    refine_repeatmasker,
)
from drayte.reporting.SummaryFilesGen import run_summary_files
from drayte.storage.archive import ArchiveStore
from drayte.storage.policy import cleanup_stage_scratch, compact_completed_extension, compact_run, output_lock, storage_mode

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="DRayTE pipeline runner")
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    parser.add_argument(
        "--repeatmodeler-library",
        default=None,
        help=(
            "Existing RepeatModeler consensi.fa.classified. Overrides the "
            "repeatmodeler_library value in the config file and skips "
            "RepeatModeler discovery."
        ),
    )
    parser.add_argument(
        "--mmseqs-rescue-db",
        default=None,
        help=(
            "Reference FASTA for classification rescue. Overrides "
            "classification.mmseqs_rescue_db in the config file."
        ),
    )
    parser.add_argument(
        "--storage-mode", choices=("keep", "compact"), default=None,
        help="Override storage.mode: retain loose evidence, or archive it after successful completion.",
    )
    return parser.parse_args()

#def parse_args() -> argparse.Namespace:
#    parser = argparse.ArgumentParser(description="DRayTE pipeline runner")
#    parser.add_argument("--config", required=True, help="Path to YAML config file")
#    return parser.parse_args()


def write_manifest(outdir: Path, name: str, data: dict) -> None:
    outpath = outdir / f"{name}.manifest.json"
    with open(outpath, "w") as handle:
        json.dump(data, handle, indent=2)


def refinement_checkpoint_is_current(
    required_outputs: list[Path],
    manifest_path: Path,
    expected_parameters: dict,
) -> bool:
    if not all(path.exists() and path.stat().st_size > 0 for path in required_outputs):
        return False
    if not manifest_path.exists() or manifest_path.stat().st_size == 0:
        return False

    try:
        cached = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False

    if cached.get("refinement_algorithm_version") != REFINEMENT_ALGORITHM_VERSION:
        return False

    return all(cached.get(key) == value for key, value in expected_parameters.items())


def log_final_results(
    logger,
    classification_result: dict,
    final_annotation_result: dict,
    refinement_result: dict,
    summary_no_filter_result: dict,
    summary_filtered_result: dict,
    family_inspection_result: dict,
) -> None:
    """Log a result map after both fresh and resumed runs."""
    logger.info("=" * 80)
    logger.info("FINAL RESULTS")

    def log_section(title: str, result: dict, fields: list[tuple[str, str]]) -> None:
        logger.info("")
        logger.info("%s", title)
        if result.get("outdir"):
            logger.info("  Output directory: %s", result["outdir"])
        for label, key in fields:
            value = result.get(key)
            if value:
                logger.info("  %s: %s", label, value)

    log_section(
        "1. CLASSIFIED TE LIBRARY AND EVIDENCE",
        classification_result,
        [
            ("Final classified library", "classified_library"),
            ("Family classifications", "classifications_tsv"),
            ("Classification evidence", "evidence_tsv"),
            ("Non-autonomous candidates", "nonauto_candidates_tsv"),
        ],
    )
    logger.info(
        "  Families: %s total; %s classified; %s unknown",
        classification_result.get("n_families", "unavailable"),
        classification_result.get("n_classified", "unavailable"),
        classification_result.get("n_unknown", "unavailable"),
    )

    log_section(
        "2. FINAL REPEATMASKER GENOME ANNOTATION",
        final_annotation_result,
        [
            ("RepeatMasker annotation", "repeatmasker_out"),
            ("RepeatMasker GFF", "repeatmasker_gff"),
            ("RepeatMasker summary table", "repeatmasker_tbl"),
            ("Masked genome", "masked_genome"),
            ("RepeatMasker alignments", "align"),
            ("RepeatMasker execution log", "repeatmasker_log"),
        ],
    )

    log_section(
        "3. REFINED AND FILTERED ANNOTATION",
        refinement_result,
        [
            ("All refined loci (TSV)", "refined_tsv"),
            ("All refined loci (GFF3)", "refined_gff"),
            ("All refined loci (BED)", "refined_bed"),
            ("Filtered non-overlapping loci (TSV)", "filtered_tsv"),
            ("Filtered non-overlapping loci (GFF3)", "filtered_gff"),
            ("Filtered non-overlapping loci (BED)", "filtered_bed"),
            ("Nested loci (GFF3)", "nested_gff"),
            ("Refinement statistics", "stats_tsv"),
        ],
    )

    log_section(
        "4. TE SUMMARY - UNFILTERED REFINED LOCI",
        summary_no_filter_result,
        [
            ("High-level count table", "high_level_count"),
            ("Family-level count table", "family_level_count"),
            ("Divergence summary", "divergence_summary"),
            ("Defragmented repeat table", "defragmented_repeats"),
        ],
    )
    logger.info("  Plots and PDF summaries: %s", summary_no_filter_result["outdir"])

    log_section(
        "5. TE SUMMARY - FILTERED NON-OVERLAPPING LOCI",
        summary_filtered_result,
        [
            ("High-level count table", "high_level_count"),
            ("Family-level count table", "family_level_count"),
            ("Divergence summary", "divergence_summary"),
            ("Defragmented repeat table", "defragmented_repeats"),
        ],
    )
    logger.info("  Plots and PDF summaries: %s", summary_filtered_result["outdir"])

    if family_inspection_result.get("enabled", True):
        logger.info("")
        logger.info("6. PER-FAMILY INSPECTION")
        logger.info("  Output directory: %s", family_inspection_result["outdir"])
        logger.info(
            "  Family manifest: %s",
            family_inspection_result.get("manifest", "unavailable"),
        )
        logger.info(
            "  Families represented: %s",
            family_inspection_result.get("families", "unavailable"),
        )
        logger.info(
            "  Per-family evidence, plots, and archives: %s",
            family_inspection_result["outdir"],
        )
    else:
        logger.info("")
        logger.info("6. PER-FAMILY INSPECTION: disabled")

    logger.info("=" * 80)


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    if args.repeatmodeler_library:
        config.extra["repeatmodeler_library"] = str(
            Path(args.repeatmodeler_library).expanduser().resolve()
        )
    if args.mmseqs_rescue_db is not None:
        classification_cfg = config.extra.get("classification")
        if classification_cfg is None:
            classification_cfg = {}
            config.extra["classification"] = classification_cfg
        if not isinstance(classification_cfg, dict):
            raise ValueError(
                "The classification config section must be a mapping "
                "to use --mmseqs-rescue-db"
            )
        classification_cfg["mmseqs_rescue_db"] = str(
            Path(args.mmseqs_rescue_db).expanduser().resolve()
        )
    if args.storage_mode is not None:
        settings = config.extra.get("storage")
        if settings is None:
            settings = {}
            config.extra["storage"] = settings
        if not isinstance(settings, dict):
            raise ValueError("The storage config section must be a mapping")
        settings["mode"] = args.storage_mode
    mode = storage_mode(config)
    # Restoration precedes every stage, including its checkpoint validation.
    with output_lock(config.outdir_path):
        ArchiveStore(config.outdir_path).restore_all()
        run_pipeline(config, args.config)
        if mode == "compact":
            compact_run(config, setup_logger(log_file=str(config.outdir_path / "drayte.log")),
                        runtime_inputs_known=True)


def run_pipeline(config, config_file=None) -> None:
    ensure_dir(config.outdir_path)
    logger = setup_logger(log_file=str(config.outdir_path / "drayte.log"))

    if config_file is not None:
        logger.info("Loaded config from %s", config_file)
    logger.info("Running pipeline for species=%s", config.species)

	# Run RepeatModeler
    discovery_result = discovery.run(config, logger)
    write_manifest(config.outdir_path / "discovery", "discovery", discovery_result)
    cleanup_stage_scratch(config, "discovery", logger)

	# Run Extention - Part 1
    extension_result = extension.run(config, discovery_result, logger)
    write_manifest(config.outdir_path / "extension", "extension", extension_result)
    compact_completed_extension(config, extension_result, logger)

	# Run Re-classify RepeatModeler libray using extentions - Part 2
    reclassify_result = reclassify.run(config, extension_result, logger)
    write_manifest(config.outdir_path / "reclassify", "reclassify", reclassify_result)

	# Run curation - Part 3
    curation_result = curation.run(config, reclassify_result, logger)
    write_manifest(config.outdir_path / "curation", "curation", curation_result)

	# Optional - Run HELIANO to identify full length Helitron 1 and 2
    heliano_result = heliano.run(config, curation_result, logger)
    write_manifest(config.outdir_path / "heliano", "heliano", heliano_result)
    cleanup_stage_scratch(config, "heliano", logger)

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
    cleanup_stage_scratch(config, "classification", logger)
# you can comment from this point on if you want to test the previous step
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
    refinement_manifest = (
        refinement_outdir / f"{config.species}.annotation_refinement.manifest.json"
    )

    refinement_cfg = config.extra.get("annotation_refinement", {})

    def refinement_option(name: str, legacy_name: str, default):
        if legacy_name in config.extra:
            return config.extra[legacy_name]
        return refinement_cfg.get(name, default)

    refinement_parameters = {
        "max_gap": int(refinement_option("max_gap", "refinement_max_gap", 150)),
        "include_nested": bool(
            refinement_option("include_nested", "refinement_include_nested", False)
        ),
        "mode": refinement_option("mode", "refinement_mode", "loose"),
        "max_consensus_overlap": int(
            refinement_option(
                "max_consensus_overlap",
                "refinement_max_consensus_overlap",
                50,
            )
        ),
        "max_consensus_jump": int(
            refinement_option(
                "max_consensus_jump",
                "refinement_max_consensus_jump",
                10000,
            )
        ),
        "resolve_overlaps": bool(
            refinement_option(
                "resolve_overlaps",
                "refinement_resolve_overlaps",
                True,
            )
        ),
        "overlap_score": refinement_option(
            "overlap_score",
            "refinement_overlap_score",
            "longest_lowdiv",
        ),
    }

    if refinement_checkpoint_is_current(
        required_outputs=[refined_tsv, filtered_tsv, stats_tsv],
        manifest_path=refinement_manifest,
        expected_parameters=refinement_parameters,
    ):
        logger.info("=" * 80)
        logger.info("STAGE: Defragmenting_RModAnnotation")
        logger.info("Output directory: %s", refinement_outdir)
        logger.info("=" * 80)
        logger.info("TE-refine outputs already exist; skipping")

        refinement_result = json.loads(refinement_manifest.read_text())
    else:
        if any(path.exists() for path in [refined_tsv, filtered_tsv, stats_tsv]):
            logger.info(
                "Existing refinement checkpoint is outdated or incompatible; regenerating"
            )
        refinement_result = refine_repeatmasker(
            rmout=Path(final_annotation_result["repeatmasker_out"]),
            species=config.species,
            outdir=refinement_outdir,
            **refinement_parameters,
            logger=logger,
        )

    refinement_result.setdefault("outdir", str(refinement_outdir))

    write_manifest(
        refinement_outdir,
        "Defragmenting_RModAnnotation",
        refinement_result,
    )

    refinement_outdir = config.outdir_path / "Defragmenting_RModAnnotation"

    summary_no_filtering_outdir = (
        config.outdir_path.parent / f"{config.species}.Defragmenting_RModAnnotation_noFiltering"
    )

    summary_no_filter_result = run_summary_files(
        refined_tsv=refinement_outdir / f"{config.species}.annotation_refinement.tsv",
        genome=Path(config.genome),
        species=config.species,
        outdir=config.outdir_path / f"{config.species}.Defragmenting_RModAnnotation_noFiltered",
        max_merge_gap=int(config.extra.get("summary_max_merge_gap", 100)),
        min_nested_overlap_fraction=float(config.extra.get("summary_min_nested_overlap_fraction", 0.80)),
        logger=logger,
    )
    
    summary_filtered_result = run_summary_files(
        refined_tsv=refinement_outdir / f"{config.species}.filteredRepeats.tsv",
        genome=Path(config.genome),
        species=config.species,
        outdir=config.outdir_path / f"{config.species}.Defragmenting_RModAnnotation_Filtered",
        max_merge_gap=int(config.extra.get("summary_max_merge_gap", 100)),
        min_nested_overlap_fraction=float(config.extra.get("summary_min_nested_overlap_fraction", 0.80)),
        logger=logger,
    )

    write_manifest(
        refinement_outdir,
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

    log_final_results(
        logger=logger,
        classification_result=classification_result,
        final_annotation_result=final_annotation_result,
        refinement_result=refinement_result,
        summary_no_filter_result=summary_no_filter_result,
        summary_filtered_result=summary_filtered_result,
        family_inspection_result=family_inspection_result,
    )
    logger.info("Pipeline completed successfully")

if __name__ == "__main__":
    main()
