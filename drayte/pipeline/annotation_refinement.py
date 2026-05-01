from __future__ import annotations

import json
from pathlib import Path

from drayte.refinement.refine_repeatmasker import refine_repeatmasker
from drayte.utils.paths import ensure_dir, stage_dir


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("ok\n")


def validate_input(path: Path, label: str) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")


def run(config, final_annotation_result: dict, logger) -> dict:
    species = config.species
    outdir = ensure_dir(stage_dir(config.outdir_path, "annotation_refinement"))
    ckpt = ensure_dir(outdir / "checkpoints")
    manifest = outdir / "annotation_refinement.manifest.json"

    max_gap = int(config.extra.get("annotation_refinement", {}).get("max_gap", 150))
    include_nested = bool(
        config.extra.get("annotation_refinement", {}).get("include_nested", False)
    )

    rmout = Path(final_annotation_result["repeatmasker_out"])

    logger.info("=" * 80)
    logger.info("STAGE: annotation_refinement")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    if (ckpt / "03_outputs.validated").exists() and manifest.exists():
        logger.info("Annotation refinement already validated; skipping")
        return json.loads(manifest.read_text())

    if not (ckpt / "01_inputs.ok").exists():
        validate_input(rmout, "RepeatMasker .out")
        touch(ckpt / "01_inputs.ok")

    result = refine_repeatmasker(
        rmout=rmout,
        species=species,
        outdir=outdir,
        max_gap=max_gap,
        include_nested=include_nested,
    )

    touch(ckpt / "02_refinement.done")
    touch(ckpt / "03_outputs.validated")

    result = {
        "stage": "annotation_refinement",
        "outdir": str(outdir),
        **result,
    }

    write_json(manifest, result)
    logger.info(
        "Annotation refinement completed: %s raw hits -> %s refined loci; %s nested hits",
        result["raw_hits"],
        result["refined_loci"],
        result["nested_hits"],
    )

    return result
