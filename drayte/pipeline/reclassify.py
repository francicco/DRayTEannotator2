from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from drayte.utils.paths import stage_dir, ensure_dir


def gather_extended_consensi(final_consensus_dir: Path, taxon: str, outfile: Path, logger) -> int:
    files = sorted(final_consensus_dir.glob(f"{taxon}*.fa"))
    if not files:
        files = sorted(final_consensus_dir.glob("*.fa"))

    if not files:
        raise FileNotFoundError(f"No consensus FASTA files found in {final_consensus_dir}")

    logger.info("Merging %d consensus FASTA files into %s", len(files), outfile)

    with open(outfile, "w") as out_handle:
        for fasta in files:
            with open(fasta) as in_handle:
                shutil.copyfileobj(in_handle, out_handle)

    return len(files)


def run_repeatclassifier(consensi_fa: Path, outdir: Path, repeatclassifier_bin: str, logger) -> Path:
    outdir = ensure_dir(outdir)
    expected = outdir / f"{consensi_fa.name}.classified"

    if expected.exists() and expected.stat().st_size > 0:
        logger.info("RepeatClassifier output already exists: %s", expected)
        return expected

    cmd = [
        repeatclassifier_bin,
        "-consensi",
        str(consensi_fa),
    ]

    logger.info("Running RepeatClassifier")
    process = subprocess.Popen(
        cmd,
        cwd=str(outdir),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert process.stdout is not None
    for line in process.stdout:
        line = line.rstrip()
        if line:
            logger.info("[reclassify] %s", line)

    rc = process.wait()
    if rc != 0:
        raise RuntimeError(f"RepeatClassifier failed with exit code {rc}")

    if not expected.exists():
        raise FileNotFoundError(
            f"RepeatClassifier finished but expected output was not found: {expected}"
        )

    return expected


def run(config, extension_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "reclassify")

    logger.info("=" * 80)
    logger.info("STAGE: reclassify")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    final_consensus_dir = Path(extension_result["extended_library_dir"])
    merged_consensi = outdir / f"{config.species}_extended_rep.fa"

    n_files = gather_extended_consensi(
        final_consensus_dir=final_consensus_dir,
        taxon=config.species,
        outfile=merged_consensi,
        logger=logger,
    )

    repeatclassifier_bin = config.extra.get("repeatclassifier_bin", "RepeatClassifier")
    classified = run_repeatclassifier(
        consensi_fa=merged_consensi,
        outdir=outdir,
        repeatclassifier_bin=repeatclassifier_bin,
        logger=logger,
    )

    result = {
        "stage": "reclassify",
        "outdir": str(outdir),
        "merged_library": str(merged_consensi),
        "classified_library": str(classified),
        "n_input_fastas": n_files,
    }

    logger.info("Reclassify stage completed")
    return result
