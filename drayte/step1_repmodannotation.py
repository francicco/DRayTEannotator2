#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from Bio import SeqIO


LOGGER = logging.getLogger("drayte.step1")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run RepeatModeler on a genome assembly, normalize family headers, "
            "and run a first RepeatMasker pass using the RepeatModeler library."
        )
    )

    parser.add_argument("--genome", required=True, help="Input genome FASTA (.fa/.fasta, optionally .gz)")
    parser.add_argument("--outdir", required=True, help="Output directory for Step1")
    parser.add_argument("--species", required=True, help="Short species prefix, e.g. Cant")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument(
        "--repeatmodeler-dir",
        required=True,
        help="Directory containing RepeatModeler executables (BuildDatabase, RepeatModeler)",
    )
    parser.add_argument(
        "--repeatscout-dir",
        required=True,
        help="RepeatScout installation root containing RepeatScout and build_lmer_table",
    )
    parser.add_argument(
        "--repeatmasker-bin",
        default="RepeatMasker",
        help="RepeatMasker executable path or name",
    )
    parser.add_argument(
        "--log-file",
        default=None,
        help="Optional log file path; default is <outdir>/step1.log",
    )

    return parser.parse_args()


def setup_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_command(cmd: list[str], cwd: Optional[Path] = None, logger: Optional[logging.Logger] = None) -> None:
    logger = logger or LOGGER
    logger.info("Running: %s", " ".join(map(str, cmd)))

    process = subprocess.Popen(
        [str(x) for x in cmd],
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    assert process.stdout is not None
    for line in process.stdout:
        line = line.rstrip()
        if line:
            logger.info(line)

    rc = process.wait()
    if rc != 0:
        raise RuntimeError(f"Command failed with exit code {rc}: {' '.join(map(str, cmd))}")


def run_command_to_logger_and_file(
    cmd: list[str],
    log_file: Path,
    cwd: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
    prefix: str = "subprocess",
) -> None:
    logger = logger or LOGGER
    logger.info("Running: %s", " ".join(map(str, cmd)))

    with open(log_file, "w") as logh:
        process = subprocess.Popen(
            [str(x) for x in cmd],
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )

        assert process.stdout is not None
        for line in process.stdout:
            line = line.rstrip()
            if not line:
                continue
            logger.info("[%s] %s", prefix, line)
            logh.write(line + "\n")
            logh.flush()

        rc = process.wait()
        if rc != 0:
            raise RuntimeError(f"Command failed with exit code {rc}: {' '.join(map(str, cmd))}")


def validate_inputs(
    genome_path: Path,
    repeatmodeler_dir: Path,
    repeatscout_dir: Path,
) -> None:
    if not genome_path.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {genome_path}")

    builddatabase_bin = repeatmodeler_dir / "BuildDatabase"
    repeatmodeler_bin = repeatmodeler_dir / "RepeatModeler"

    if not builddatabase_bin.exists():
        raise FileNotFoundError(f"BuildDatabase not found: {builddatabase_bin}")
    if not repeatmodeler_bin.exists():
        raise FileNotFoundError(f"RepeatModeler not found: {repeatmodeler_bin}")

    expected_rs = [
        repeatscout_dir / "RepeatScout",
        repeatscout_dir / "build_lmer_table",
    ]
    missing = [str(p) for p in expected_rs if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Invalid RepeatScout directory; missing expected files: " + ", ".join(missing)
        )


def prepare_genome(genome_path: Path, assemblies_dir: Path, species: str, logger: Optional[logging.Logger] = None) -> Path:
    logger = logger or LOGGER
    dest = assemblies_dir / f"{species}.fa"

    if dest.exists() and dest.stat().st_size > 0:
        logger.info("Genome FASTA already exists: %s", dest)
        return dest

    logger.info("Preparing genome FASTA: %s -> %s", genome_path, dest)

    if genome_path.suffix == ".gz":
        with gzip.open(genome_path, "rt") as fin, open(dest, "wt") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copy2(genome_path, dest)

    return dest


def repeatmodeler_db_exists(assemblies_dir: Path, species: str) -> bool:
    candidates = [
        assemblies_dir / f"{species}.nhr",
        assemblies_dir / f"{species}.nin",
        assemblies_dir / f"{species}.nsq",
    ]
    return all(p.exists() for p in candidates)


def find_repeatmodeler_library(rmodeler_dir: Path) -> Optional[Path]:
    local_lib = rmodeler_dir / "consensi.fa.classified"
    if local_lib.exists() and local_lib.stat().st_size > 0:
        return local_lib

    rm_dirs = sorted(
        [p for p in rmodeler_dir.iterdir() if p.is_dir() and p.name.startswith("RM_")],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )

    for rm_dir in rm_dirs:
        candidate = rm_dir / "consensi.fa.classified"
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate

    return None


def parse_repeatmodeler_output_dir(rm_log: Path) -> Optional[Path]:
    if not rm_log.exists():
        return None

    with open(rm_log) as handle:
        for line in handle:
            if line.startswith("Using output directory"):
                try:
                    return Path(line.strip().split("=", 1)[1].strip())
                except IndexError:
                    return None
    return None


def log_tail(log_file: Path, n: int = 80, logger: Optional[logging.Logger] = None) -> None:
    logger = logger or LOGGER
    if not log_file.exists():
        logger.error("Log file does not exist: %s", log_file)
        return

    try:
        with open(log_file) as handle:
            lines = handle.readlines()[-n:]
        logger.error("Last %d lines of %s:", n, log_file)
        for line in lines:
            logger.error(line.rstrip())
    except Exception as exc:
        logger.error("Could not read log tail from %s: %s", log_file, exc)


def copy_if_missing(src: Path, dst: Path, logger: Optional[logging.Logger] = None) -> None:
    logger = logger or LOGGER
    if dst.exists() and dst.stat().st_size > 0:
        logger.info("File already exists: %s", dst)
        return
    shutil.copyfile(src, dst)


def normalize_family_headers(
    original_fasta: Path,
    edited_fasta: Path,
    species: str,
    logger: Optional[logging.Logger] = None,
) -> None:
    logger = logger or LOGGER
    logger.info("Normalizing RepeatModeler family headers: %s -> %s", original_fasta, edited_fasta)

    with open(original_fasta) as infile, open(edited_fasta, "w") as outfile:
        for i, record in enumerate(SeqIO.parse(infile, "fasta"), start=1):
            raw_id = record.id.split()[0]

            if "#" in raw_id:
                name, cls = raw_id.split("#", 1)
                name = name.replace("rnd-", f"{species}-rnd-")
                new_id = f"{name}#{cls}"
            else:
                name = raw_id.replace("rnd-", f"{species}-rnd-")
                new_id = name

            record.id = new_id
            record.description = new_id
            SeqIO.write(record, outfile, "fasta")


def repeatmasker_outputs_exist(assemblies_dir: Path, rmasker_dir: Path, species: str) -> bool:
    candidates = [
        assemblies_dir / f"{species}.fa.tbl",
        rmasker_dir / f"{species}.fa.tbl",
    ]
    return any(p.exists() and p.stat().st_size > 0 for p in candidates)


def move_repeatmasker_outputs(
    assemblies_dir: Path,
    rmasker_dir: Path,
    species: str,
    logger: Optional[logging.Logger] = None,
) -> None:
    logger = logger or LOGGER
    for suffix in [".out", ".tbl", ".cat.gz", ".masked", ".gff"]:
        src = assemblies_dir / f"{species}.fa{suffix}"
        dst = rmasker_dir / src.name
        if src.exists():
            if dst.exists():
                logger.info("RepeatMasker output already moved: %s", dst)
            else:
                logger.info("Moving %s -> %s", src, dst)
                shutil.move(str(src), str(dst))


def run_step1(
    genome: Path,
    outdir: Path,
    species: str,
    threads: int,
    repeatmodeler_dir: Path,
    repeatscout_dir: Path,
    repeatmasker_bin: str = "RepeatMasker",
    logger: Optional[logging.Logger] = None,
) -> dict:
    logger = logger or LOGGER

    outdir = outdir.resolve()
    genome = genome.resolve()
    repeatmodeler_dir = repeatmodeler_dir.resolve()
    repeatscout_dir = repeatscout_dir.resolve()

    assemblies_dir = ensure_dir(outdir / "assemblies_dir")
    rmodeler_dir = ensure_dir(outdir / "rmodeler_dir")
    rmasker_dir = ensure_dir(outdir / "rmasker_dir")

    validate_inputs(genome, repeatmodeler_dir, repeatscout_dir)

    builddatabase_bin = repeatmodeler_dir / "BuildDatabase"
    repeatmodeler_bin = repeatmodeler_dir / "RepeatModeler"

    logger.info("Starting Step1.RepModAnnotation")
    logger.info("Genome: %s", genome)
    logger.info("Species: %s", species)
    logger.info("Threads: %d", threads)
    logger.info("Output directory: %s", outdir)

    logger.info("Step 1/4: preparing genome FASTA")
    genome_fa = prepare_genome(genome, assemblies_dir, species, logger=logger)

    logger.info("Step 2/4: checking/building RepeatModeler database")
    db_prefix = assemblies_dir / species
    rm_log = rmodeler_dir / f"{species}.RMrun.out"

    if repeatmodeler_db_exists(assemblies_dir, species):
        logger.info("RepeatModeler database already exists for %s", db_prefix)
    else:
        run_command(
            [
                str(builddatabase_bin),
                "-name", str(db_prefix),
                str(genome_fa),
            ],
            logger=logger,
        )

    logger.info("Step 3/4: checking/running RepeatModeler")
    consensi = find_repeatmodeler_library(rmodeler_dir)

    if consensi is None:
        logger.info("No existing RepeatModeler library found; running RepeatModeler")
        try:
            run_command_to_logger_and_file(
                [
                    str(repeatmodeler_bin),
                    "-rscout_dir", str(repeatscout_dir),
                    "-database", str(db_prefix),
                    "-threads", str(threads),
                ],
                log_file=rm_log,
                cwd=rmodeler_dir,
                logger=logger,
                prefix="RepeatModeler",
            )
        except RuntimeError:
            log_tail(rm_log, n=80, logger=logger)
            raise RuntimeError(f"RepeatModeler failed; see log: {rm_log}")

        consensi = find_repeatmodeler_library(rmodeler_dir)

        if consensi is None:
            rmoutdir = parse_repeatmodeler_output_dir(rm_log)
            if rmoutdir is not None:
                candidate = rmoutdir / "consensi.fa.classified"
                if candidate.exists() and candidate.stat().st_size > 0:
                    consensi = candidate

        if consensi is None:
            raise RuntimeError(
                "RepeatModeler appears to have completed, but no valid consensi.fa.classified was found."
            )
    else:
        logger.info("Found existing RepeatModeler library: %s", consensi)

    logger.info("Using RepeatModeler library: %s", consensi)

    copied_consensi = rmodeler_dir / "consensi.fa.classified"
    if consensi.resolve() != copied_consensi.resolve():
        copy_if_missing(consensi, copied_consensi, logger=logger)
    else:
        logger.info("Library already in standard location: %s", copied_consensi)

    species_families = rmodeler_dir / f"{species}-families.fa"
    copy_if_missing(copied_consensi, species_families, logger=logger)

    edited_fasta = rmodeler_dir / f"{species}-families.mod.fa"
    if edited_fasta.exists() and edited_fasta.stat().st_size > 0:
        logger.info("Normalized family FASTA already exists: %s", edited_fasta)
    else:
        normalize_family_headers(copied_consensi, edited_fasta, species, logger=logger)

    logger.info("Step 4/4: checking/running RepeatMasker")
    if repeatmasker_outputs_exist(assemblies_dir, rmasker_dir, species):
        logger.info("RepeatMasker outputs already exist; skipping RepeatMasker")
    else:
        run_command(
            [
                str(repeatmasker_bin),
                "-xsmall",
                "-gc",
                "-gff",
                "-norna",
                "-pa", str(threads),
                "-lib", str(copied_consensi),
                "-s",
                str(genome_fa),
            ],
            logger=logger,
        )

    move_repeatmasker_outputs(assemblies_dir, rmasker_dir, species, logger=logger)

    logger.info("Step1 completed successfully")
    logger.info("Key outputs:")
    logger.info("  Raw library: %s", copied_consensi)
    logger.info("  Species family FASTA: %s", species_families)
    logger.info("  Normalized library: %s", edited_fasta)
    logger.info("  RepeatMasker dir: %s", rmasker_dir)

    return {
        "outdir": str(outdir),
        "assemblies_dir": str(assemblies_dir),
        "rmodeler_dir": str(rmodeler_dir),
        "rmasker_dir": str(rmasker_dir),
        "genome_fa": str(genome_fa),
        "raw_library": str(copied_consensi),
        "species_families": str(species_families),
        "normalized_library": str(edited_fasta),
        "repeatmasker_tbl": str(rmasker_dir / f"{species}.fa.tbl"),
    }


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir).resolve()
    log_file = Path(args.log_file).resolve() if args.log_file else outdir / "step1.log"
    setup_logging(log_file)

    run_step1(
        genome=Path(args.genome).resolve(),
        outdir=outdir,
        species=args.species,
        threads=args.threads,
        repeatmodeler_dir=Path(args.repeatmodeler_dir).resolve(),
        repeatscout_dir=Path(args.repeatscout_dir).resolve(),
        repeatmasker_bin=args.repeatmasker_bin,
        logger=LOGGER,
    )


if __name__ == "__main__":
    main()
