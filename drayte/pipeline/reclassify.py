from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO

from drayte.utils.paths import stage_dir, ensure_dir


REPEATCLASSIFIER_CHECKPOINT_SCHEMA_VERSION = 1


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


def fasta_sequences_by_family(fasta: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for record in SeqIO.parse(str(fasta), "fasta"):
        family_id = record.id.split("#", 1)[0]
        if family_id in sequences:
            raise ValueError(f"Duplicate family identifier in {fasta}: {family_id}")
        sequences[family_id] = str(record.seq).upper()
    return sequences


def compare_classified_output(
    consensi_fa: Path,
    classified_fa: Path,
) -> dict[str, int] | None:
    if not classified_fa.exists() or classified_fa.stat().st_size == 0:
        return None
    try:
        input_sequences = fasta_sequences_by_family(consensi_fa)
        output_sequences = fasta_sequences_by_family(classified_fa)
    except (OSError, ValueError):
        return None
    input_ids = set(input_sequences)
    output_ids = set(output_sequences)
    shared = input_ids & output_ids
    return {
        "input_families": len(input_ids),
        "output_families": len(output_ids),
        "missing_families": len(input_ids - output_ids),
        "unexpected_families": len(output_ids - input_ids),
        "changed_sequences": sum(
            input_sequences[family_id] != output_sequences[family_id]
            for family_id in shared
        ),
    }


def classified_output_matches_input(consensi_fa: Path, classified_fa: Path) -> bool:
    comparison = compare_classified_output(consensi_fa, classified_fa)
    return bool(
        comparison
        and comparison["input_families"] > 0
        and comparison["missing_families"] == 0
        and comparison["unexpected_families"] == 0
        and comparison["changed_sequences"] == 0
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repeatclassifier_checkpoint_path(outdir: Path) -> Path:
    return outdir / "repeatclassifier.checkpoint.json"


def write_json_atomic(path: Path, data: dict) -> None:
    temporary = Path(f"{path}.partial")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def checkpoint_data(
    consensi_fa: Path,
    classified_fa: Path,
    repeatclassifier_bin: str,
) -> dict:
    return {
        "schema_version": REPEATCLASSIFIER_CHECKPOINT_SCHEMA_VERSION,
        "status": "complete",
        "input_fasta": str(consensi_fa.resolve()),
        "input_sha256": sha256_file(consensi_fa),
        "output_fasta": str(classified_fa.resolve()),
        "output_sha256": sha256_file(classified_fa),
        "family_count": len(fasta_sequences_by_family(consensi_fa)),
        "repeatclassifier_bin": str(repeatclassifier_bin),
    }


def checkpoint_is_valid(
    checkpoint: Path,
    consensi_fa: Path,
    classified_fa: Path,
    repeatclassifier_bin: str,
) -> bool:
    try:
        cached = json.loads(checkpoint.read_text())
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        return False
    if not isinstance(cached, dict):
        return False
    try:
        current = checkpoint_data(
            consensi_fa,
            classified_fa,
            repeatclassifier_bin,
        )
    except (FileNotFoundError, OSError, ValueError):
        return False
    return cached == current


def legacy_manifest_proves_completion(
    outdir: Path,
    consensi_fa: Path,
    classified_fa: Path,
) -> bool:
    manifest = outdir / "reclassify.manifest.json"
    try:
        cached = json.loads(manifest.read_text())
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        return False
    if not isinstance(cached, dict) or cached.get("stage") != "reclassify":
        return False
    try:
        cached_input = Path(cached["merged_library"]).resolve()
        cached_output = Path(cached["classified_library"]).resolve()
    except (KeyError, TypeError):
        return False
    return (
        cached_input == consensi_fa.resolve()
        and cached_output == classified_fa.resolve()
        and manifest.stat().st_mtime_ns >= classified_fa.stat().st_mtime_ns
    )


def archive_stale_output(path: Path) -> Path:
    candidate = path.with_name(f"{path.name}.stale")
    counter = 1
    while candidate.exists():
        candidate = path.with_name(f"{path.name}.stale.{counter}")
        counter += 1
    path.replace(candidate)
    return candidate


def run_repeatclassifier(consensi_fa: Path, outdir: Path, repeatclassifier_bin: str, logger) -> Path:
    outdir = ensure_dir(outdir)
    expected = outdir / f"{consensi_fa.name}.classified"
    checkpoint = repeatclassifier_checkpoint_path(outdir)

    if classified_output_matches_input(consensi_fa, expected):
        if checkpoint_is_valid(
            checkpoint,
            consensi_fa,
            expected,
            repeatclassifier_bin,
        ):
            logger.info(
                "Validated completed RepeatClassifier checkpoint: %s",
                expected,
            )
            return expected
        if legacy_manifest_proves_completion(outdir, consensi_fa, expected):
            write_json_atomic(
                checkpoint,
                checkpoint_data(consensi_fa, expected, repeatclassifier_bin),
            )
            logger.info(
                "Validated and migrated completed legacy RepeatClassifier "
                "checkpoint: %s",
                expected,
            )
            return expected
        archived = archive_stale_output(expected)
        logger.warning(
            "Existing RepeatClassifier output contains the current family "
            "sequences but has no successful completion checkpoint; moved "
            "it to %s",
            archived,
        )

    elif expected.exists():
        comparison = compare_classified_output(consensi_fa, expected)
        archived = archive_stale_output(expected)
        if comparison is None:
            logger.warning(
                "Existing RepeatClassifier output is empty or invalid; "
                "moved it to %s",
                archived,
            )
        else:
            logger.warning(
                "Existing RepeatClassifier output does not match the current "
                "merged library (%d missing families, %d unexpected families, "
                "%d changed sequences); moved it to %s",
                comparison["missing_families"],
                comparison["unexpected_families"],
                comparison["changed_sequences"],
                archived,
            )

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

    if not classified_output_matches_input(consensi_fa, expected):
        comparison = compare_classified_output(consensi_fa, expected)
        detail = "missing or invalid output"
        if comparison is not None:
            detail = (
                f"{comparison['missing_families']} missing families, "
                f"{comparison['unexpected_families']} unexpected families, "
                f"{comparison['changed_sequences']} changed sequences"
            )
        raise RuntimeError(
            "RepeatClassifier exited successfully but its output does not "
            f"match the merged input ({detail}): {expected}"
        )

    write_json_atomic(
        checkpoint,
        checkpoint_data(consensi_fa, expected, repeatclassifier_bin),
    )
    logger.info("Recorded completed RepeatClassifier checkpoint: %s", checkpoint)
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
