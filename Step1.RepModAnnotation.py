#!/usr/bin/env python3

import argparse
import gzip
import logging
import shutil
import subprocess
from pathlib import Path
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
        help="Reserved for compatibility; currently not used directly.",
    )
    parser.add_argument(
        "--repeatmodeler-dir",
        required=True,
        help="Directory containing RepeatModeler executables (BuildDatabase, RepeatModeler).",
    )
    parser.add_argument(
        "--repeatscout-dir",
        required=True,
        help="Directory containing RepeatScout binaries.",
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

    handlers = [
        logging.FileHandler(log_file),
        logging.StreamHandler(),
    ]

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        handlers=handlers,
    )


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_command(cmd, cwd: Path | None = None) -> None:
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


def prepare_genome(genome_path: Path, assemblies_dir: Path, species: str) -> Path:
    """
    Copy or decompress the input genome into assemblies_dir as <species>.fa
    """
    dest = assemblies_dir / f"{species}.fa"

    if dest.exists():
        LOGGER.info("Genome FASTA already exists: %s", dest)
        return dest

    LOGGER.info("Preparing genome FASTA: %s -> %s", genome_path, dest)

    if genome_path.suffix == ".gz":
        with gzip.open(genome_path, "rt") as fin, open(dest, "wt") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copy2(genome_path, dest)

    return dest


def parse_repeatmodeler_output_dir(rm_log: Path) -> Path:
    with open(rm_log, "r") as handle:
        for line in handle:
            if line.startswith("Using output directory"):
                return Path(line.strip().split("=")[1].strip())

    raise RuntimeError(f"Could not find RepeatModeler output directory in {rm_log}")


def normalize_family_headers(original_fasta: Path, edited_fasta: Path, species: str) -> None:
    LOGGER.info("Normalizing RepeatModeler family headers")
    with open(original_fasta) as infile, open(edited_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            header = record.id.split()[0]
            header = header.replace("#", "__")
            header = header.replace("/", "___")
            header = header.replace("rnd", f"{species}-rnd")
            record.id = header
            record.description = header
            SeqIO.write(record, outfile, "fasta")


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir).resolve()
    assemblies_dir = ensure_dir(outdir / "assemblies_dir")
    rmodeler_dir = ensure_dir(outdir / "rmodeler_dir")
    rmasker_dir = ensure_dir(outdir / "rmasker_dir")

    log_file = Path(args.log_file) if args.log_file else outdir / "step1.log"
    setup_logging(log_file)

    genome_path = Path(args.genome).resolve()
    repeatmodeler_dir = Path(args.repeatmodeler_dir).resolve()
    repeatscout_dir = Path(args.repeatscout_dir).resolve()

    builddatabase_bin = repeatmodeler_dir / "BuildDatabase"
    repeatmodeler_bin = repeatmodeler_dir / "RepeatModeler"

    LOGGER.info("Starting Step1.RepModAnnotation")
    LOGGER.info("Genome: %s", genome_path)
    LOGGER.info("Species: %s", args.species)
    LOGGER.info("Threads: %s", args.threads)
    LOGGER.info("Output directory: %s", outdir)

    genome_fa = prepare_genome(genome_path, assemblies_dir, args.species)

    db_prefix = assemblies_dir / args.species
    rm_log = rmodeler_dir / f"{args.species}.RMrun.out"

    # Build RepeatModeler database
    if not (assemblies_dir / f"{args.species}.nhr").exists() and not (assemblies_dir / f"{args.species}.nin").exists():
        run_command([
            builddatabase_bin,
            "-name", db_prefix,
            genome_fa,
        ])
    else:
        LOGGER.info("RepeatModeler database appears to exist already for %s", db_prefix)

    # Run RepeatModeler
    if not rm_log.exists():
        LOGGER.info("Running RepeatModeler")
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
            raise RuntimeError(f"RepeatModeler failed; see log: {rm_log}")
    else:
        LOGGER.info("RepeatModeler log exists already: %s", rm_log)

    rmoutdir = parse_repeatmodeler_output_dir(rm_log)
    LOGGER.info("RepeatModeler output directory: %s", rmoutdir)

    consensi = rmoutdir / "consensi.fa.classified"
    if not consensi.exists():
        raise FileNotFoundError(f"Expected RepeatModeler library not found: {consensi}")

    copied_consensi = rmodeler_dir / "consensi.fa.classified"
    if not copied_consensi.exists():
        shutil.copyfile(consensi, copied_consensi)

    species_families = rmodeler_dir / f"{args.species}-families.fa"
    if not species_families.exists():
        shutil.copyfile(consensi, species_families)

    edited_fasta = rmodeler_dir / f"{args.species}-families.mod.fa"
    if not edited_fasta.exists():
        normalize_family_headers(consensi, edited_fasta, args.species)
    else:
        LOGGER.info("Normalized family FASTA exists already: %s", edited_fasta)

    # First-pass RepeatMasker
    rm_tbl = assemblies_dir / f"{args.species}.fa.tbl"
    if not rm_tbl.exists():
        run_command([
            args.repeatmasker_bin,
            "-xsmall",
            "-gc",
            "-gff",
            "-norna",
            "-pa", str(args.threads),
            "-lib", str(consensi),
            "-s",
            str(genome_fa),
        ])
    else:
        LOGGER.info("RepeatMasker output already exists: %s", rm_tbl)

    # Move RepeatMasker outputs into rmasker_dir if present
    for suffix in [".out", ".tbl", ".cat.gz", ".masked", ".gff"]:
        src = assemblies_dir / f"{args.species}.fa{suffix}"
        if src.exists():
            dst = rmasker_dir / src.name
            if not dst.exists():
                shutil.move(str(src), str(dst))

    LOGGER.info("Step1 completed successfully")
    LOGGER.info("Key outputs:")
    LOGGER.info("  Raw library: %s", copied_consensi)
    LOGGER.info("  Normalized library: %s", edited_fasta)
    LOGGER.info("  RepeatMasker dir: %s", rmasker_dir)


if __name__ == "__main__":
    main()
