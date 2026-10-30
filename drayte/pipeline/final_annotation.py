from __future__ import annotations

import shutil
from pathlib import Path

from drayte.utils.paths import stage_dir, ensure_dir
from drayte.utils.subprocess import run_command_to_file
from drayte.utils.repeatmasker_checkpoint import (
    begin_attempt, checkpoint_is_valid, checkpoint_path, input_signature,
    legacy_outputs_usable, remove_stale_source_outputs, validate_outputs,
    write_checkpoint,
)


def final_annotation_exists(outdir: Path, taxon: str) -> bool:
    required = [
        outdir / f"{taxon}.fa.out",
        outdir / f"{taxon}.fa.tbl",
        outdir / f"{taxon}.fa.masked",
        outdir / f"{taxon}.fa.cat.gz",
        outdir / f"{taxon}.fa.align",
    ]

    return all(f.exists() and f.stat().st_size > 0 for f in required)


def move_repeatmasker_outputs(genome_fa: Path, outdir: Path, logger) -> None:
    parent = genome_fa.parent
    genome_name = genome_fa.name

    for suffix in [".out", ".tbl", ".gff", ".masked", ".cat.gz", ".align"]:
        src = parent / f"{genome_name}{suffix}"
        dst = outdir / f"{genome_name}{suffix}"

        if src.exists():
            if dst.exists():
                dst.unlink()
            logger.info("Moving %s -> %s", src, dst)
            shutil.move(str(src), str(dst))


def run(config, curation_result: dict, logger) -> dict:
    outdir = ensure_dir(stage_dir(config.outdir_path, "final_annotation"))

    logger.info("=" * 80)
    logger.info("STAGE: final_annotation")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    taxon = config.species
    genome_fa = config.outdir_path / "discovery" / "assemblies_dir" / f"{taxon}.fa"

    if "classified_library" not in curation_result:
        raise FileNotFoundError(
            "classified_library missing from final_annotation input. "
            "Run classification before final_annotation."
        )
    
    final_library = Path(curation_result["classified_library"])
    
    if not final_library.exists() or final_library.stat().st_size == 0:
        raise FileNotFoundError(f"Classified library not found: {final_library}")    
    
    tbl = outdir / f"{taxon}.fa.tbl"
    out = outdir / f"{taxon}.fa.out"
    gff = outdir / f"{taxon}.fa.gff"
    masked = outdir / f"{taxon}.fa.masked"
    cat_gz = outdir / f"{taxon}.fa.cat.gz"
    align = outdir / f"{taxon}.fa.align"
    repeatmasker_log = outdir / f"{taxon}.RepeatMasker.log"

    outputs_exist = final_annotation_exists(outdir, taxon)
    required = [out, tbl, masked, cat_gz, align]
    checkpoint = checkpoint_path(outdir)
    options = [str(config.extra.get("repeatmasker_bin", "RepeatMasker")),
               "-s", "-xsmall", "-gc", "-gff", "-norna", "-a"]
    signature = input_signature(genome_fa, final_library, options)
    outputs_current = checkpoint_is_valid(checkpoint, signature, required)
    if (not checkpoint.exists() and outputs_exist
            and out.stat().st_mtime >= final_library.stat().st_mtime
            and legacy_outputs_usable(required, repeatmasker_log)):
        write_checkpoint(checkpoint, signature, required, legacy=True)
        outputs_current = True
        logger.warning("Adopted complete legacy RepeatMasker outputs; historical input provenance is unavailable")
    
    if outputs_current:
        logger.info("Validated completed final RepeatMasker checkpoint; skipping")
    else:
        if outputs_exist:
            logger.info(
                "Final RepeatMasker checkpoint or inputs changed; re-running"
            )
        cmd = [
            str(config.extra.get("repeatmasker_bin", "RepeatMasker")),
            "-s",
            "-xsmall",
            "-gc",
            "-gff",
            "-norna",
            "-a",
            "-pa", str(config.threads),
            "-lib", str(final_library),
            str(genome_fa),
        ]

        begin_attempt(checkpoint, signature)
        remove_stale_source_outputs(genome_fa, [".out", ".tbl", ".gff", ".masked", ".cat.gz", ".align"])
        run_command_to_file(
            cmd,
            log_file=repeatmasker_log,
            logger=logger,
            cwd=genome_fa.parent,
            label="RepeatMasker",
        )
        validate_outputs([genome_fa.parent / f"{genome_fa.name}{suffix}"
                          for suffix in (".out", ".tbl", ".masked", ".cat.gz", ".align")])
        move_repeatmasker_outputs(genome_fa, outdir, logger)
        write_checkpoint(checkpoint, signature, required)

    result = {
        "stage": "final_annotation",
        "outdir": str(outdir),
        "repeatmasker_tbl": str(tbl),
        "repeatmasker_out": str(out),
        "repeatmasker_gff": str(gff),
        "masked_genome": str(masked),
        "cat_gz": str(cat_gz),
        "align": str(align),
    }
    if repeatmasker_log.is_file():
        result["repeatmasker_log"] = str(repeatmasker_log.resolve())
    return result
