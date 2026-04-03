from __future__ import annotations

import csv
import gzip
import logging
import shutil
import subprocess
from pathlib import Path

from drayte.utils.paths import stage_dir, ensure_dir
from drayte.extension import (
    filter_fasta_by_length,
    rename_repeatmodeler_headers,
    run_extract_align,
    run_extend_consensus,
    postprocess_extension_outputs,
    categorize_extension,
)


def gunzip_to_file(src_gz: Path, dest_fa: Path, logger) -> None:
    if dest_fa.exists() and dest_fa.stat().st_size > 0:
        logger.info("Genome FASTA already present: %s", dest_fa)
        return

    logger.info("Decompressing genome: %s -> %s", src_gz, dest_fa)
    with gzip.open(src_gz, "rt") as fin, open(dest_fa, "wt") as fout:
        shutil.copyfileobj(fin, fout)


def run_blast(query_fasta: Path, genome_fa: Path, blast_dir: Path, threads: int, logger) -> Path:
    db_check = blast_dir / f"{genome_fa.name}.nsq"
    blast_out = blast_dir / f"{genome_fa.stem}_blastn.out"

    if not db_check.exists():
        logger.info("Building BLAST database")
        subprocess.run(
            ["makeblastdb", "-in", str(genome_fa), "-dbtype", "nucl"],
            cwd=str(blast_dir),
            check=True,
        )

    if not blast_out.exists() or blast_out.stat().st_size == 0:
        logger.info("Running BLASTN")
        subprocess.run(
            [
                "blastn",
                "-query", str(query_fasta),
                "-db", str(genome_fa),
                "-outfmt", "6",
                "-out", str(blast_out),
                "-num_threads", str(threads),
            ],
            cwd=str(blast_dir),
            check=True,
        )
    else:
        logger.info("BLAST output already exists: %s", blast_out)

    return blast_out


def ensure_twobit(genome_fa: Path, genome_2bit: Path, logger) -> None:
    if genome_2bit.exists() and genome_2bit.stat().st_size > 0:
        logger.info("2bit genome already exists: %s", genome_2bit)
        return

    logger.info("Creating 2bit genome: %s", genome_2bit)
    subprocess.run(
        ["faToTwoBit", str(genome_fa), str(genome_2bit)],
        check=True,
    )


def write_summary(summary_file: Path, rows: list[dict]) -> None:
    fieldnames = [
        "te_id",
        "cat_file",
        "rep_fa",
        "msa_fa",
        "png_file",
        "hit_count",
        "consensus_length",
        "category",
    ]
    with open(summary_file, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run(config, discovery_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "extension")

    logger.info("=" * 80)
    logger.info("STAGE: extension")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    taxon = config.species
    workdir = config.outdir_path
    assemblies_dir = workdir / "discovery" / "assemblies_dir"
    rmodeler_dir = workdir / "discovery" / "rmodeler_dir"

    extensions_dir = ensure_dir(outdir)
    extensionwork = ensure_dir(extensions_dir / "extensionwork")
    extendlogs = ensure_dir(extensions_dir / "extendlogs")
    blastfiles = ensure_dir(extensions_dir / "blastfiles")
    extract_dir = ensure_dir(extensions_dir / "extract_align")
    images = ensure_dir(extensions_dir / "images_and_alignments")
    rejects = ensure_dir(images / "rejects")
    likely_tes = ensure_dir(images / "likely_TEs")
    possible_sd = ensure_dir(images / "possible_SD")
    final_consensuses = ensure_dir(extensions_dir / "final_consensuses")

    genome_fa = assemblies_dir / f"{taxon}.fa"
    genome_2bit = assemblies_dir / f"{taxon}.2bit"

    genome_input = config.genome_path.resolve()
    if genome_input.suffix == ".gz":
        gunzip_to_file(genome_input, genome_fa, logger)

    ensure_twobit(genome_fa, genome_2bit, logger)

    raw_consensus = Path(discovery_result["raw_library"])
    filtered_fasta = rmodeler_dir / f"{taxon}-families.filtered.fa"
    final_query_fasta = rmodeler_dir / f"{taxon}-families.mod.fa"

    if not filtered_fasta.exists():
        n_kept = filter_fasta_by_length(raw_consensus, filtered_fasta, min_length=100)
        logger.info("Filtered consensuses >=100 bp: %d", n_kept)

    if not final_query_fasta.exists():
        n_renamed = rename_repeatmodeler_headers(filtered_fasta, final_query_fasta, species=taxon)
        logger.info("Renamed consensus headers: %d", n_renamed)

    blast_query_copy = blastfiles / final_query_fasta.name
    if not blast_query_copy.exists():
        shutil.copy2(final_query_fasta, blast_query_copy)

    blast_out = run_blast(
        query_fasta=blast_query_copy,
        genome_fa=genome_fa,
        blast_dir=blastfiles,
        threads=config.threads,
        logger=logger,
    )

    extracted = run_extract_align(
        genome_fasta=genome_fa,
        blast_file=blast_out,
        library_fasta=final_query_fasta,
        output_dir=extract_dir,
        max_hits_per_query=int(config.extra.get("extension_max_hits", 50)),
        flank_left=int(config.extra.get("extension_left_buffer", 100)),
        flank_right=int(config.extra.get("extension_right_buffer", 100)),
    )

    extend_script = config.extra["repeatmodeler_extend_script"]
    summary_rows = []

    for te_id, cat_file in extracted.items():
        rep_file = final_consensuses / f"{te_id}_rep.fa"
        te_workdir = extensionwork / te_id
        te_log = extendlogs / f"{te_id}.extend.log"

        if not rep_file.exists():
            run_extend_consensus(
                extend_script=extend_script,
                genome_2bit=genome_2bit,
                family_fasta=cat_file,
                workdir=te_workdir,
                logfile=te_log,
            )

        rep_fa, msa_fa, png_file, hit_count, consensus_length = postprocess_extension_outputs(
            te_id=te_id,
            te_workdir=te_workdir,
        )

        category = categorize_extension(hit_count, consensus_length)

        if rep_fa and not rep_file.exists():
            shutil.copy2(rep_fa, rep_file)

        if category == "reject":
            if msa_fa:
                shutil.copy2(msa_fa, rejects / msa_fa.name)
            if png_file:
                shutil.copy2(png_file, rejects / png_file.name)
        elif category == "possible_SD":
            if rep_fa:
                shutil.copy2(rep_fa, possible_sd / rep_fa.name)
            if msa_fa:
                shutil.copy2(msa_fa, possible_sd / msa_fa.name)
            if png_file:
                shutil.copy2(png_file, possible_sd / png_file.name)
        else:
            if rep_fa:
                shutil.copy2(rep_fa, likely_tes / rep_fa.name)
            if msa_fa:
                shutil.copy2(msa_fa, likely_tes / msa_fa.name)
            if png_file:
                shutil.copy2(png_file, likely_tes / png_file.name)

        summary_rows.append(
            {
                "te_id": te_id,
                "cat_file": str(cat_file),
                "rep_fa": str(rep_fa) if rep_fa else "",
                "msa_fa": str(msa_fa) if msa_fa else "",
                "png_file": str(png_file) if png_file else "",
                "hit_count": hit_count,
                "consensus_length": consensus_length,
                "category": category,
            }
        )

    summary_file = extensions_dir / "extension_summary.tsv"
    write_summary(summary_file, summary_rows)

    result = {
        "stage": "extension",
        "outdir": str(extensions_dir),
        "extended_library_dir": str(final_consensuses),
        "summary_tsv": str(summary_file),
    }

    logger.info("Extension stage completed")
    return result
