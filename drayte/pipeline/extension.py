from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from drayte.utils.paths import stage_dir, ensure_dir
from drayte.extension import (
    filter_fasta_by_length,
    rename_repeatmodeler_headers,
    run_extract_align,
    run_extend_consensus,
    postprocess_extension_outputs,
    categorize_extension,
)

from drayte.extension.extract_align import safe_filename

def gunzip_to_file(src_gz: Path, dest_fa: Path, logger) -> None:
    if dest_fa.exists() and dest_fa.stat().st_size > 0:
        logger.info("Genome FASTA already present: %s", dest_fa)
        return

    logger.info("Decompressing genome: %s -> %s", src_gz, dest_fa)
    with gzip.open(src_gz, "rt") as fin, open(dest_fa, "wt") as fout:
        shutil.copyfileobj(fin, fout)


def ensure_twobit(genome_fa: Path, genome_2bit: Path, logger) -> None:
    if genome_2bit.exists() and genome_2bit.stat().st_size > 0:
        logger.info("2bit genome already exists: %s", genome_2bit)
        return

    logger.info("Creating 2bit genome: %s", genome_2bit)
    subprocess.run(
        ["faToTwoBit", str(genome_fa), str(genome_2bit)],
        check=True,
    )


def merge_candidate_libraries(
    primary_fasta: Path,
    secondary_fasta: Optional[Path],
    merged_fasta: Path,
    logger,
) -> Path:

    if merged_fasta.exists() and merged_fasta.stat().st_size > 0:
        logger.info("Merged candidate library already exists: %s", merged_fasta)
        return merged_fasta

    seen = set()
    records = []

    for fasta in [primary_fasta, secondary_fasta]:
        if fasta is None:
            continue
        if not fasta.exists() or fasta.stat().st_size == 0:
            continue

        for rec in SeqIO.parse(str(fasta), "fasta"):
            if rec.id in seen:
                continue
            seen.add(rec.id)
            records.append(rec)

    with open(merged_fasta, "w") as handle:
        SeqIO.write(records, handle, "fasta")

    logger.info("Merged candidate libraries into %s (%d records)", merged_fasta, len(records))
    return merged_fasta


def run_blast(query_fasta: Path, genome_fa: Path, threads: int, logger) -> Path:
    """
    BLAST DB is built next to genome_fa because makeblastdb uses the genome path as DB prefix.
    BLAST output is stored next to query_fasta.
    """
    db_files = [Path(str(genome_fa) + ext) for ext in [".nhr", ".nin", ".nsq"]]
    blast_out = query_fasta.parent / f"{genome_fa.stem}_blastn.out"

    if not all(p.exists() and p.stat().st_size > 0 for p in db_files):
        logger.info("Building BLAST database for %s", genome_fa)
        subprocess.run(
            ["makeblastdb", "-in", str(genome_fa), "-dbtype", "nucl"],
            check=True,
        )
    else:
        logger.info("BLAST database already exists for %s", genome_fa)

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
            check=True,
        )
    else:
        logger.info("BLAST output already exists: %s", blast_out)

    return blast_out


def write_summary(summary_file: Path, rows: list[dict]) -> None:
    fieldnames = [
        "te_id",
        "cat_file",
        "status",
        "category",
        "rep_fa",
        "msa_fa",
        "png_file",
        "hit_count",
        "consensus_length",
        "error_message",
    ]
    with open(summary_file, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


#def run(config, discovery_result: dict, structure_result: Optional[dict], logger) -> dict:
def run(config, discovery_result: dict, logger) -> dict:
    outdir = stage_dir(config.outdir_path, "extension")

    logger.info("=" * 80)
    logger.info("STAGE: extension")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    taxon = config.species
    discovery_dir = config.outdir_path / "discovery"
    assemblies_dir = discovery_dir / "assemblies_dir"
    rmodeler_dir = discovery_dir / "rmodeler_dir"

    extensionwork = ensure_dir(outdir / "extensionwork")
    extendlogs = ensure_dir(outdir / "extendlogs")
    blastfiles = ensure_dir(outdir / "blastfiles")
    extract_dir = ensure_dir(outdir / "extract_align")
    images = ensure_dir(outdir / "images_and_alignments")
    rejects = ensure_dir(images / "rejects")
    likely_tes = ensure_dir(images / "likely_TEs")
    possible_sd = ensure_dir(images / "possible_SD")
    final_consensuses = ensure_dir(outdir / "final_consensuses")

    genome_fa = assemblies_dir / f"{taxon}.fa"
    genome_2bit = assemblies_dir / f"{taxon}.2bit"

    genome_input = config.genome_path.resolve()
    if genome_input.suffix == ".gz":
        gunzip_to_file(genome_input, genome_fa, logger)

    ensure_twobit(genome_fa, genome_2bit, logger)

    raw_consensus = Path(discovery_result["raw_library"])
    candidate_library = raw_consensus

    filtered_fasta = rmodeler_dir / f"{taxon}-families.filtered.fa"
    final_query_fasta = rmodeler_dir / f"{taxon}-families.mod.fa"

    if not filtered_fasta.exists() or filtered_fasta.stat().st_size == 0:
        n_kept = filter_fasta_by_length(candidate_library, filtered_fasta, min_length=100)
        logger.info("Filtered consensuses >=100 bp: %d", n_kept)
    else:
        logger.info("Filtered candidate FASTA already exists: %s", filtered_fasta)

    if not final_query_fasta.exists() or final_query_fasta.stat().st_size == 0:
        n_renamed = rename_repeatmodeler_headers(filtered_fasta, final_query_fasta, species=taxon)
        logger.info("Renamed consensus headers: %d", n_renamed)
    else:
        logger.info("Normalized candidate FASTA already exists: %s", final_query_fasta)

    blast_query_copy = blastfiles / final_query_fasta.name
    if not blast_query_copy.exists() or blast_query_copy.stat().st_size == 0:
        shutil.copy2(final_query_fasta, blast_query_copy)

    blast_out = run_blast(
        query_fasta=blast_query_copy,
        genome_fa=genome_fa,
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
    summary_rows: list[dict] = []

    for te_id, cat_file in extracted.items():
        safe_te_id = safe_filename(te_id)
        
        rep_file = final_consensuses / f"{safe_te_id}_rep.fa"
        te_workdir = extensionwork / safe_te_id
        te_log = extendlogs / f"{safe_te_id}.extend.log"

        rep_fa = None
        msa_fa = None
        png_file = None
        hit_count = 0
        consensus_length = 0
        category = "unknown"
        status = "ok"
        error_message = ""

        try:
            if not rep_file.exists() or rep_file.stat().st_size == 0:
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

            if not rep_fa or not Path(rep_fa).exists() or Path(rep_fa).stat().st_size == 0:
                status = "failed_extension"
                error_message = "No rep.fa produced"
                raise RuntimeError(error_message)

            if rep_fa and rep_fa.exists():
                if not rep_file.exists() or rep_file.stat().st_size == 0:
                    shutil.copy2(rep_fa, rep_file)

            target_dir = (
                rejects if category == "reject"
                else possible_sd if category == "possible_SD"
                else likely_tes
            )

            if rep_fa and rep_fa.exists():
                dst = target_dir / rep_fa.name
                if not dst.exists():
                    shutil.copy2(rep_fa, dst)

            if msa_fa and msa_fa.exists():
                dst = target_dir / msa_fa.name
                if not dst.exists():
                    shutil.copy2(msa_fa, dst)

            if png_file and png_file.exists():
                dst = target_dir / png_file.name
                if not dst.exists():
                    shutil.copy2(png_file, dst)

        except Exception as exc:
            status = "failed_extension"
            error_message = str(exc)
            logger.warning("Could not extend %s: %s", cat_file, exc)

        summary_rows.append(
            {
                "te_id": te_id,
                "cat_file": str(cat_file),
                "status": status,
                "category": category,
                "rep_fa": str(rep_fa) if rep_fa else "",
                "msa_fa": str(msa_fa) if msa_fa else "",
                "png_file": str(png_file) if png_file else "",
                "hit_count": hit_count,
                "consensus_length": consensus_length,
                "error_message": error_message,
            }
        )

    summary_file = outdir / "extension_summary.tsv"
    write_summary(summary_file, summary_rows)

    n_ok = sum(1 for row in summary_rows if row["status"] == "ok")
    n_fail = sum(1 for row in summary_rows if row["status"] != "ok")
    logger.info("Extension stage completed: %d ok, %d failed", n_ok, n_fail)

    result = {
        "stage": "extension",
        "outdir": str(outdir),
        "extended_library_dir": str(final_consensuses),
        "summary_tsv": str(summary_file),
        "n_ok": n_ok,
        "n_failed": n_fail,
    }

    return result

