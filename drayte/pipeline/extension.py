from __future__ import annotations

import json
import csv
import gzip
import hashlib
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def extension_input_signature(
    query_fasta: Path,
    genome_fasta: Path,
    config,
) -> dict:
    extend_script = Path(config.extra["repeatmodeler_extend_script"]).resolve()
    return {
        "schema_version": 2,
        "species": config.species,
        "query_fasta_sha256": sha256_file(query_fasta),
        "genome_fasta_sha256": sha256_file(genome_fasta),
        "extension_max_hits": int(config.extra.get("extension_max_hits", 50)),
        "extension_left_buffer": int(config.extra.get("extension_left_buffer", 100)),
        "extension_right_buffer": int(config.extra.get("extension_right_buffer", 100)),
        "repeatmodeler_extend_script": str(extend_script),
        "repeatmodeler_extend_script_sha256": sha256_file(extend_script),
        "ucsc_tools_dir": str(Path(config.extra["ucsc_tools_dir"]).resolve()),
        "repeatmodeler_dir": str(Path(config.extra["repeatmodeler_dir"]).resolve()),
        "repeat_after_me_dir": str(
            Path(config.extra["repeat_after_me_dir"]).resolve()
        ),
        "phrap_dir": str(Path(config.extra["phrap_dir"]).resolve()),
    }


def extension_signatures_compatible(cached: dict, current: dict) -> bool:
    if cached.get("schema_version") == current.get("schema_version"):
        return cached == current

    if cached.get("schema_version") == 1 and current.get("schema_version") == 2:
        version_one_keys = {
            "schema_version",
            "species",
            "query_fasta_sha256",
            "genome_fasta_sha256",
            "extension_max_hits",
            "extension_left_buffer",
            "extension_right_buffer",
            "repeatmodeler_extend_script",
        }
        normalized_current = {
            key: (1 if key == "schema_version" else current.get(key))
            for key in version_one_keys
        }
        return all(cached.get(key) == value for key, value in normalized_current.items())

    return False


def validate_extension_outputs(
    result: dict,
    extensionwork: Path,
    final_consensuses: Path,
) -> tuple[bool, str]:
    summary_value = result.get("summary_tsv")
    if not summary_value:
        return False, "manifest has no summary_tsv"

    summary_file = Path(summary_value)
    if not summary_file.exists() or summary_file.stat().st_size == 0:
        return False, f"missing or empty summary: {summary_file}"

    with open(summary_file) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    n_ok = sum(row.get("status") == "ok" for row in rows)
    n_failed = len(rows) - n_ok
    if n_ok != int(result.get("n_ok", -1)):
        return False, "summary n_ok does not match manifest"
    if n_failed != int(result.get("n_failed", -1)):
        return False, "summary n_failed does not match manifest"
    if not rows:
        return False, "summary contains no families"
    if n_failed:
        return False, f"summary contains {n_failed} failed families"

    for row in rows:
        if row.get("status") != "ok":
            continue

        te_id = row.get("te_id", "")
        safe_te_id = safe_filename(te_id)
        workdir = extensionwork / safe_te_id
        final_rep = final_consensuses / f"{safe_te_id}_rep.fa"

        if not workdir.is_dir():
            return False, f"missing extension work directory: {workdir}"
        if not final_rep.exists() or final_rep.stat().st_size == 0:
            return False, f"missing or empty final consensus: {final_rep}"

    return True, "outputs are complete"


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

    if not final_query_fasta.exists() or final_query_fasta.stat().st_size == 0:
        if not filtered_fasta.exists() or filtered_fasta.stat().st_size == 0:
            n_kept = filter_fasta_by_length(
                candidate_library,
                filtered_fasta,
                min_length=100,
            )
            logger.info("Filtered consensuses >=100 bp: %d", n_kept)

        n_renamed = rename_repeatmodeler_headers(
            filtered_fasta,
            final_query_fasta,
            species=taxon,
        )
        logger.info("Renamed consensus headers: %d", n_renamed)

    input_signature = extension_input_signature(
        query_fasta=final_query_fasta,
        genome_fasta=genome_fa,
        config=config,
    )

    manifest = outdir / "extension.manifest.json"
    if manifest.exists():
        try:
            cached_result = json.loads(manifest.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning("Could not read extension manifest; resuming stage: %s", exc)
        else:
            cached_signature = cached_result.get("input_signature")
            if cached_signature is not None and not extension_signatures_compatible(
                cached_signature,
                input_signature,
            ):
                raise RuntimeError(
                    "Existing extension outputs were created from different "
                    "inputs or settings. Use a new outdir, or move the existing "
                    "extension and downstream stage directories aside before "
                    "restarting."
                )

            blast_query_copy = blastfiles / final_query_fasta.name
            if cached_signature is None:
                if (
                    not blast_query_copy.exists()
                    or blast_query_copy.stat().st_size == 0
                    or sha256_file(blast_query_copy)
                    != input_signature["query_fasta_sha256"]
                ):
                    raise RuntimeError(
                        "Legacy extension manifest cannot be matched to the "
                        "current query library. Use a new outdir, or move the "
                        "existing extension and downstream stage directories "
                        "aside before restarting."
                    )

            valid, reason = validate_extension_outputs(
                cached_result,
                extensionwork=extensionwork,
                final_consensuses=final_consensuses,
            )
            if valid:
                cached_result["input_signature"] = input_signature
                logger.info("Validated existing extension outputs; skipping")
                return cached_result

            logger.warning("Extension cache is incomplete (%s); resuming", reason)

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
    if blast_query_copy.exists() and blast_query_copy.stat().st_size > 0:
        if sha256_file(blast_query_copy) != input_signature["query_fasta_sha256"]:
            raise RuntimeError(
                "Cached extension query library differs from the current "
                "library. Use a new outdir, or move the existing extension "
                "and downstream stage directories aside before restarting."
            )
    else:
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

    ucsc_tools_dir = config.extra["ucsc_tools_dir"]
    repeatmodeler_dir = config.extra["repeatmodeler_dir"]
    repeat_after_me_dir = config.extra["repeat_after_me_dir"]
    phrap_dir = config.extra["phrap_dir"]

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
                    ucsc_tools_dir=ucsc_tools_dir,
                    repeatmodeler_dir=repeatmodeler_dir,
                    repeat_after_me_dir=repeat_after_me_dir,
                    phrap_dir=phrap_dir,
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
        "input_signature": input_signature,
    }

    if n_ok == 0:
        manifest.write_text(json.dumps(result, indent=2) + "\n")
        raise RuntimeError(
            "Extension produced no successful consensuses "
            f"(0 ok, {n_fail} failed); see {summary_file} and the extendlogs directory"
        )

    return result

