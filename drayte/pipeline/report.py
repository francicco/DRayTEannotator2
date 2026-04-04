from __future__ import annotations

from pathlib import Path

from drayte.utils.paths import stage_dir, ensure_dir


def run(config, curation_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "report")
    ensure_dir(outdir)

    logger.info("Starting report stage")

    report_txt = outdir / "final_outputs.txt"
    with open(report_txt, "w") as handle:
        handle.write(f"curated_library\t{curation_result.get('curated_library', '')}\n")
        handle.write(f"curated_full\t{curation_result.get('curated_full', '')}\n")
        handle.write(f"curated_nr\t{curation_result.get('curated_nr', '')}\n")
        handle.write(f"curated_metadata\t{curation_result.get('curated_metadata', '')}\n")
        handle.write(f"run_summary\t{curation_result.get('run_summary', '')}\n")

    result = {
        "stage": "report",
        "outdir": str(outdir),
        "report_txt": str(report_txt),
        "curated_library": curation_result.get("curated_library", ""),
    }

    logger.info("Report stage completed")
    return result
