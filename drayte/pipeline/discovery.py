from drayte.utils.paths import stage_dir
from drayte.utils.subprocess import run_command


def run(config, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "discovery")

    logger.info("=" * 80)
    logger.info("STAGE: discovery")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    cmd = [
        "python",
        config.tools.step1_repmodannotation,
        "--genome", str(config.genome_path.resolve()),
        "--outdir", str(outdir),
        "--species", config.species,
        "--threads", str(config.threads),
        "--batches", str(config.extra.get("batches", 1)),
        "--repeatmodeler-dir", config.extra["repeatmodeler_dir"],
        "--repeatscout-dir", config.extra["repeatscout_dir"],
        "--repeatmasker-bin", config.extra.get("repeatmasker_bin", "RepeatMasker"),
    ]

    run_command(cmd, logger, prefix="discovery")

    result = {
        "stage": "discovery",
        "outdir": str(outdir),
        "genome": str(config.genome_path.resolve()),
        "raw_library": str(outdir / "rmodeler_dir" / "consensi.fa.classified"),
        "normalized_library": str(outdir / "rmodeler_dir" / f"{config.species}-families.mod.fa"),
        "repeatmasker_dir": str(outdir / "rmasker_dir"),
    }

    logger.info("Discovery stage completed")
    return result
