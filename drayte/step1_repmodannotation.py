#!/usr/bin/env python3

import argparse
import gzip
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from Bio import SeqIO


LOGGER = logging.getLogger("Step1.RepModAnnotation")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run RepeatModeler on a genome assembly, normalize family headers, "
            "and run a first RepeatMasker pass using the RepeatModeler library."
        )
    )

    parser.add_argument(
        "--genome",
        required=True,
        help="Path to input genome FASTA (.fa, .fasta, optionally gzipped).",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Working/output directory for this step.",
    )
    parser.add_argument(
        "--species",
        required=True,
        help="Short species/genome prefix to use in output file names.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads to use.",
    )
    parser.add_argument(
        "--batches",
        type=int,
        default=1,
        help="Reserved for compatibility; currently unused.",
    )
    parser.add_argument(
        "--repeatmodeler-dir",
        required=True,
        help="Directory containing RepeatModeler executables (BuildDatabase, RepeatModeler).",
    )
    parser.add_argument(
        "--repeatscout-dir",
        required=True,
        help="Directory containing RepeatScout installation root.",
    )
    parser.add_argument(
        "--repeatmasker-bin",
        default="RepeatMasker",
        help="RepeatMasker executable name or full path.",
    )
    parser.add_argument(
        "--log-file",
        default=None,
        help="Optional log file path. Defaults to <outdir>/step1.log",
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


def run_command(cmd: list[str], cwd: Optional[Path] = None) -> None:
    LOGGER.info("Running: %s", " ".join(map(str, cmd)))
    result = subprocess.run(
        [str(x) for x in cmd],
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
    )

    if result.stdout:
        LOGGER.info(result.stdout.strip())

    if result.returncode != 0:
        if result.stderr:
            LOGGER.error(result.stderr.strip())
        raise RuntimeError(f"Command failed: {' '.join(map(str, cmd))}")


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
            "Invalid RepeatScout directory; missing expected files: "
            + ", ".join(missing)
        )


def prepare_genome(genome_path: Path, assemblies_dir: Path, species: str) -> Path:
    """
    Copy or decompress the input genome into assemblies_dir as <species>.fa
    """
    dest = assemblies_dir / f"{species}.fa"

    if dest.exists() and dest.stat().st_size > 0:
        LOGGER.info("Genome FASTA already exists: %s", dest)
        return dest

    LOGGER.info("Preparing genome FASTA: %s -> %s", genome_path, dest)

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
    """
    Return the best existing RepeatModeler library if present.

    Preference:
    1. rmodeler_dir/consensi.fa.classified
    2. newest RM_*/consensi.fa.classified
    """
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

    with open(rm_log, "r") as handle:
        for line in handle:
            if line.startswith("Using output directory"):
                try:
                    return Path(line.strip().split("=", 1)[1].strip())
                except IndexError:
                    return None
    return None


def log_tail(log_file: Path, n: int = 80) -> None:
    if not log_file.exists():
        LOGGER.error("Log file does not exist: %s", log_file)
        return

    try:
        with open(log_file, "r") as handle:
            lines = handle.readlines()[-n:]
        LOGGER.error("Last %d lines of %s:", n, log_file)
        for line in lines:
            LOGGER.error(line.rstrip())
    except Exception as exc:
        LOGGER.error("Could not read log tail from %s: %s", log_file, exc)


def copy_if_missing(src: Path, dst: Path) -> None:
    if dst.exists() and dst.stat().st_size > 0:
        LOGGER.info("File already exists: %s", dst)
        return
    shutil.copyfile(src, dst)


def normalize_family_headers(original_fasta: Path, edited_fasta: Path, species: str) -> None:
    LOGGER.info("Normalizing RepeatModeler family headers: %s -> %s", original_fasta, edited_fasta)
    with open(original_fasta) as infile, open(edited_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            header = record.id.split()[0]
            header = header.replace("#", "__")
            header = header.replace("/", "___")
            header = header.replace("rnd", f"{species}-rnd")
            record.id = header
            record.description = header
            SeqIO.write(record, outfile, "fasta")


def repeatmasker_outputs_exist(assemblies_dir: Path, rmasker_dir: Path, species: str) -> bool:
    targets = [
        assemblies_dir / f"{species}.fa.tbl",
        rmasker_dir / f"{species}.fa.tbl",
    ]
    return any(p.exists() and p.stat().st_size > 0 for p in targets)


def move_repeatmasker_outputs(assemblies_dir: Path, rmasker_dir: Path, species: str) -> None:
    for suffix in [".out", ".tbl", ".cat.gz", ".masked", ".gff"]:
        src = assemblies_dir / f"{species}.fa{suffix}"
        dst = rmasker_dir / src.name
        if src.exists():
            if dst.exists():
                LOGGER.info("RepeatMasker output already moved: %s", dst)
            else:
                LOGGER.info("Moving %s -> %s", src, dst)
                shutil.move(str(src), str(dst))


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir).resolve()
    assemblies_dir = ensure_dir(outdir / "assemblies_dir")
    rmodeler_dir = ensure_dir(outdir / "rmodeler_dir")
    rmasker_dir = ensure_dir(outdir / "rmasker_dir")

    log_file = Path(args.log_file).resolve() if args.log_file else outdir / "step1.log"
    setup_logging(log_file)

    genome_path = Path(args.genome).resolve()
    repeatmodeler_dir = Path(args.repeatmodeler_dir).resolve()
    repeatscout_dir = Path(args.repeatscout_dir).resolve()

    validate_inputs(genome_path, repeatmodeler_dir, repeatscout_dir)

    builddatabase_bin = repeatmodeler_dir / "BuildDatabase"
    repeatmodeler_bin = repeatmodeler_dir / "RepeatModeler"

    LOGGER.info("Starting Step1.RepModAnnotation")
    LOGGER.info("Genome: %s", genome_path)
    LOGGER.info("Species: %s", args.species)
    LOGGER.info("Threads: %d", args.threads)
    LOGGER.info("Output directory: %s", outdir)
    LOGGER.info("Step 1/4: preparing genome FASTA")
    
    genome_fa = prepare_genome(genome_path, assemblies_dir, args.species)
    db_prefix = assemblies_dir / args.species
    rm_log = rmodeler_dir / f"{args.species}.RMrun.out"

    # Build RepeatModeler database
    if repeatmodeler_db_exists(assemblies_dir, args.species):
        LOGGER.info("RepeatModeler database already exists for %s", db_prefix)
    else:
        run_command([
            str(builddatabase_bin),
            "-name", str(db_prefix),
            str(genome_fa),
        ])

    # Find or run RepeatModeler
    consensi = find_repeatmodeler_library(rmodeler_dir)

    if consensi is None:
        LOGGER.info("Step 3/4: checking/running RepeatModeler")
        LOGGER.info("No existing RepeatModeler library found; running RepeatModeler")
        with open(rm_log, "w") as log_handle:
            result = subprocess.run(
                [
                    str(repeatmodeler_bin),
                    "-rscout_dir", str(repeatscout_dir),
                    "-database", str(db_prefix),
                    "-threads", str(args.threads),
                ],
                cwd=str(rmodeler_dir),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
            )

        if result.returncode != 0:
            log_tail(rm_log, n=80)
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
        LOGGER.info("Found existing RepeatModeler library: %s", consensi)

    LOGGER.info("Using RepeatModeler library: %s", consensi)

    # Standardize key output copies in rmodeler_dir
    copied_consensi = rmodeler_dir / "consensi.fa.classified"
    if consensi.resolve() != copied_consensi.resolve():
        copy_if_missing(consensi, copied_consensi)
    else:
        LOGGER.info("Library already in standard location: %s", copied_consensi)

    species_families = rmodeler_dir / f"{args.species}-families.fa"
    copy_if_missing(copied_consensi, species_families)

    edited_fasta = rmodeler_dir / f"{args.species}-families.mod.fa"
    if edited_fasta.exists() and edited_fasta.stat().st_size > 0:
        LOGGER.info("Normalized family FASTA already exists: %s", edited_fasta)
    else:
        normalize_family_headers(copied_consensi, edited_fasta, args.species)

    # First-pass RepeatMasker
    LOGGER.info("Step 4/4: checking/running RepeatMasker")
    if repeatmasker_outputs_exist(assemblies_dir, rmasker_dir, args.species):
        LOGGER.info("RepeatMasker outputs already exist; skipping RepeatMasker")
    else:
        run_command([
            str(args.repeatmasker_bin),
            "-xsmall",
            "-gc",
            "-gff",
            "-norna",
            "-pa", str(args.threads),
            "-lib", str(copied_consensi),
            "-s",
            str(genome_fa),
        ])

    move_repeatmasker_outputs(assemblies_dir, rmasker_dir, args.species)

    LOGGER.info("Step1 completed successfully")
    LOGGER.info("Key outputs:")
    LOGGER.info("  Raw library: %s", copied_consensi)
    LOGGER.info("  Species family FASTA: %s", species_families)
    LOGGER.info("  Normalized library: %s", edited_fasta)
    LOGGER.info("  RepeatMasker dir: %s", rmasker_dir)


if __name__ == "__main__":
    main()
