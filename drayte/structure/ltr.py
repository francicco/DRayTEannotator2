from __future__ import annotations

import re
import subprocess
from pathlib import Path

from .models import StructureCandidate
from .utils import extract_region_fasta


def run_suffixerator(genome_fasta: Path, index_prefix: Path, gt_bin: str, logger) -> None:
    esq = Path(str(index_prefix) + ".esq")
    if esq.exists():
        logger.info("GenomeTools suffixerator index already exists: %s", esq)
        return

    logger.info("Building GenomeTools suffixerator index")
    subprocess.run(
        [
            gt_bin,
            "suffixerator",
            "-db", str(genome_fasta),
            "-indexname", str(index_prefix),
            "-tis", "-suf", "-lcp", "-des", "-ssp", "-sds", "-dna",
        ],
        check=True,
    )


def run_ltrharvest(index_prefix: Path, out_gff3: Path, out_fasta: Path, gt_bin: str, logger) -> None:
    if out_gff3.exists() and out_gff3.stat().st_size > 0:
        logger.info("LTRharvest output already exists: %s", out_gff3)
        return

    logger.info("Running LTRharvest")
    with open(out_fasta, "w") as fasta_handle:
        subprocess.run(
            [
                gt_bin,
                "ltrharvest",
                "-index", str(index_prefix),
                "-gff3", str(out_gff3),
                "-out", str(out_fasta),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )


def parse_ltrharvest_gff3(
    gff3_file: Path,
    genome_fasta: Path,
    candidates_dir: Path,
    species: str,
) -> list[StructureCandidate]:
    candidates: list[StructureCandidate] = []
    counter = 0

    with open(gff3_file) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue

            contig, source, feature, start, end, score, strand, phase, attrs = fields
            if feature != "LTR_retrotransposon":
                continue

            counter += 1
            start_i = int(start)
            end_i = int(end)
            candidate_id = f"{species}_LTRstruct_{counter}"
            out_fasta = candidates_dir / f"{candidate_id}.fa"

            extract_region_fasta(
                genome_fasta=genome_fasta,
                contig=contig,
                start_1based=start_i,
                end_1based=end_i,
                strand=strand,
                out_fasta=out_fasta,
                seq_id=candidate_id,
            )

            candidates.append(
                StructureCandidate(
                    candidate_id=candidate_id,
                    module="LTR",
                    contig=contig,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    length=end_i - start_i + 1,
                    fasta_path=out_fasta,
                )
            )

    return candidates


def run_ltr_module(
    genome_fasta: Path,
    outdir: Path,
    species: str,
    gt_bin: str,
    logger,
) -> list[StructureCandidate]:
    outdir.mkdir(parents=True, exist_ok=True)
    index_prefix = outdir / f"{species}.ltrindex"
    gff3_file = outdir / "ltrharvest.gff3"
    raw_fasta = outdir / "ltrharvest.raw.fa"
    candidates_dir = outdir / "candidates"
    candidates_dir.mkdir(parents=True, exist_ok=True)

    run_suffixerator(genome_fasta, index_prefix, gt_bin, logger)
    run_ltrharvest(index_prefix, gff3_file, raw_fasta, gt_bin, logger)

    candidates = parse_ltrharvest_gff3(
        gff3_file=gff3_file,
        genome_fasta=genome_fasta,
        candidates_dir=candidates_dir,
        species=species,
    )

    logger.info("LTR module recovered %d structural candidates", len(candidates))
    return candidates
