from __future__ import annotations

from pathlib import Path

from drayte.utils.paths import stage_dir, ensure_dir
from drayte.curation.inputs import prepare_curation_inputs
from drayte.curation.orfs import run_orf_detection
from drayte.curation.homology import run_repeatpep_search
from drayte.curation.tables import build_family_table
from drayte.curation.selection import select_representatives


def run(config, reclassify_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "curation")

    logger.info("=" * 80)
    logger.info("STAGE: curation")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    prepared = prepare_curation_inputs(
        classified_library=Path(reclassify_result["classified_library"]),
        outdir=outdir,
        species=config.species,
        logger=logger,
    )

    orf_result = run_orf_detection(
        input_fasta=prepared["clean_library"],
        outdir=outdir / "orfs",
        getorf_bin=config.extra.get("getorf_bin", "getorf"),
        logger=logger,
    )

    homology_result = run_repeatpep_search(
        input_fasta=prepared["clean_library"],
        outdir=outdir / "homology",
        diamond_bin=config.extra.get("diamond_bin", "diamond"),
        repeatpeps_db=config.extra.get("repeatpeps_db"),
        threads=config.threads,
        logger=logger,
    )

    family_table = build_family_table(
        library_fasta=prepared["clean_library"],
        orf_result=orf_result,
        homology_result=homology_result,
        outdir=outdir / "tables",
        logger=logger,
    )

    selected = select_representatives(
        family_table_tsv=family_table,
        library_fasta=prepared["clean_library"],
        outdir=outdir / "selected",
        logger=logger,
    )

    result = {
        "stage": "curation",
        "outdir": str(outdir),
        "family_table": str(family_table),
        "curated_library": str(selected["curated_library"]),
        "selected_table": str(selected["selected_table"]),
    }

    logger.info("Curation stage completed")
    return result
