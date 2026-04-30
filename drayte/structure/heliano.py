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

def normalize_heliano_fasta_headers(in_fasta: Path, out_fasta: Path, species: str) -> None:
    seen = {}

    with in_fasta.open() as inp, out_fasta.open("w") as out:
        for line in inp:
            if line.startswith(">"):
                raw = line[1:].strip().split()[0]

                if "HLE1_orfonly" in raw:
                    subtype = "HLE1"
                elif "HLE2_orfonly" in raw:
                    subtype = "HLE2"
                elif "HLE1_or_HLE2" in raw:
                    subtype = "HLE1_or_HLE2"
                else:
                    subtype = "Helitron"

                key = f"{species}_HELIANO_{raw}"
                seen[key] = seen.get(key, 0) + 1

                if seen[key] > 1:
                    key = f"{key}_dup{seen[key]}"

                out.write(f">{key}#RC/Helitron\n")
            else:
                out.write(line)

def heliano_outputs_exist(outdir: Path, species: str) -> bool:
    final_lib = outdir / f"{species}.heliano.fa"
    return final_lib.exists() and final_lib.stat().st_size > 0

def filter_heliano_against_curated(
    heliano_fa: Path,
    curated_fa: Path,
    outdir: Path,
    species: str,
    heliano_cfg: dict,
    logger,
) -> Path:
    import subprocess
    import shutil

    mmseqs_bin = heliano_cfg.get("mmseqs_bin", "mmseqs")
    identity = str(heliano_cfg.get("redundancy_identity", 0.90))
    query_cov = str(heliano_cfg.get("query_cov", 0.80))
    threads = str(heliano_cfg.get("threads", 1))

    tsv = outdir / f"{species}.heliano_vs_curated.tsv"
    redundant_ids = outdir / f"{species}.heliano.redundant.ids"
    unique_fa = outdir / f"{species}.heliano.unique_vs_curated.fa"
    tmpdir = outdir / "mmseqs_tmp"

    # Clean tmpdir if exists (MMseqs can fail otherwise)
    if tmpdir.exists():
        shutil.rmtree(tmpdir)

    logger.info("Running MMseqs2 redundancy filtering")

    cmd = [
        str(mmseqs_bin), "easy-search",
        str(heliano_fa),
        str(curated_fa),
        str(tsv),
        str(tmpdir),
        "--search-type", "3",            # nucleotide vs nucleotide
        "--min-seq-id", identity,
        "-c", query_cov,
        "--cov-mode", "0",               # query coverage
        "--threads", threads,
    ]

    logger.info("MMseqs command: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)

    # Parse results
    ids = set()
    if tsv.exists():
        with tsv.open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                qid = line.split("\t")[0]
                ids.add(qid)

    # Write redundant IDs
    redundant_ids.write_text("\n".join(sorted(ids)) + ("\n" if ids else ""))

    # Filter FASTA
    with heliano_fa.open() as inp, unique_fa.open("w") as out:
        write = False
        for line in inp:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                write = seq_id not in ids
            if write:
                out.write(line)

    n_unique = sum(1 for line in unique_fa.open() if line.startswith(">"))

    logger.info(
        "HELIANO redundancy filter: %s redundant, %s unique",
        len(ids),
        n_unique,
    )

    return unique_fa 

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
        logger.info("HELIANO outputs already exist; skipping HELIANO run")
    
        curated_library = Path(curation_result["final_library"])
    
        unique_lib = filter_heliano_against_curated(
            heliano_fa=final_lib,
            curated_fa=curated_library,
            outdir=outdir,
            species=species,
            heliano_cfg=heliano_cfg,
            logger=logger,
        )
    
        return {
            "stage": "heliano",
            "enabled": True,
            "outdir": str(outdir),
            "heliano_library": str(final_lib),
            "heliano_unique_library": str(unique_lib),
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
        normalize_heliano_fasta_headers(detected, final_lib, species)

    logger.info("HELIANO library: %s", final_lib)

    curated_library = Path(curation_result["final_library"])
    
    unique_lib = filter_heliano_against_curated(
        heliano_fa=final_lib,
        curated_fa=curated_library,
        outdir=outdir,
        species=species,
        heliano_cfg=heliano_cfg,
        logger=logger,
    )
    
    return {
        "stage": "heliano",
        "enabled": True,
        "outdir": str(outdir),
        "heliano_library": str(final_lib),
        "heliano_unique_library": str(unique_lib),
    }
