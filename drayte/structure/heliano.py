from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from datetime import datetime

from drayte.utils.paths import ensure_dir, stage_dir


def find_heliano_fasta(outdir: Path) -> Path | None:
    candidates = list(outdir.rglob("*.fa")) + list(outdir.rglob("*.fasta"))
    candidates = [p for p in candidates if p.exists() and p.stat().st_size > 0]

    if not candidates:
        return None

    preferred = [
        p for p in candidates
        if any(
            k in p.name.lower()
            for k in ["rep", "representative", "consensus", "helitron", "heliano"]
        )
    ]

    return preferred[0] if preferred else candidates[0]


def heliano_outputs_exist(outdir: Path, species: str) -> bool:
    final_lib = outdir / f"{species}.heliano.fa"
    return final_lib.exists() and final_lib.stat().st_size > 0


def run(config, curation_result: dict, logger) -> dict:
    species = config.species
    outdir = ensure_dir(stage_dir(config.outdir_path, "heliano"))

    heliano_cfg = config.extra.get("structure", {}).get("heliano", {})

    logger.info("=" * 80)
    logger.info("STAGE: heliano")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    if not heliano_cfg.get("enabled", False):
        logger.info("HELIANO disabled; skipping")
        return {
            "stage": "heliano",
            "enabled": False,
            "outdir": str(outdir),
            "heliano_library": None,
        }

    final_lib = outdir / f"{species}.heliano.fa"

    if heliano_outputs_exist(outdir, species):
        logger.info("HELIANO outputs already exist; skipping")
        return {
            "stage": "heliano",
            "enabled": True,
            "outdir": str(outdir),
            "heliano_library": str(final_lib),
        }

    genome_fa = config.outdir_path / "discovery" / "assemblies_dir" / f"{species}.fa"
    heliano_bin = heliano_cfg.get("bin", "heliano")

    if not genome_fa.exists() or genome_fa.stat().st_size == 0:
        raise FileNotFoundError(f"Genome FASTA not found for HELIANO: {genome_fa}")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = outdir / f"HEL_{timestamp}"

    if run_dir.exists():
        shutil.rmtree(run_dir)

    cmd = [
        str(heliano_bin),
        "-g", str(genome_fa),
        "-o", str(run_dir),
        "-n", str(config.threads),
        "--nearest",
        "-dn", "6000",
        "-flank_sim", "0.5",
        "-w", "10000",
    ]

    logger.info("Running HELIANO from cwd=%s", outdir)
    logger.info("Running HELIANO: %s", " ".join(cmd))
    subprocess.run(
        cmd,
        check=True,
        cwd=str(outdir),
    )

    detected = find_heliano_fasta(run_dir)

    if detected is None:
        logger.warning("HELIANO completed but no FASTA output was found")
        return {
            "stage": "heliano",
            "enabled": True,
            "outdir": str(outdir),
            "heliano_library": None,
        }

    if detected.resolve() != final_lib.resolve():
        shutil.copyfile(detected, final_lib)

    logger.info("HELIANO library: %s", final_lib)

    return {
        "stage": "heliano",
        "enabled": True,
        "outdir": str(outdir),
        "heliano_library": str(final_lib),
    }
