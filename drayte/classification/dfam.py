from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys

from dataclasses import dataclass
from pathlib import Path

from typing import List

from .ids import clean_family_id
from drayte.utils.fasta_split import split_fasta_balanced_by_length
from concurrent.futures import ThreadPoolExecutor, as_completed


DFAM_CHUNK_CACHE_SCHEMA_VERSION = 1
DFAM_SCAN_CACHE_SCHEMA_VERSION = 1
CHUNK_FILE_PATTERN = re.compile(r"^chunk_(\d+)\.(?:fa|tblout)$")

@dataclass
class DfamHit:
    family_id: str
    model_name: str
    accession: str
    evalue: float
    score: float
    ali_start: int
    ali_end: int
    dfam_class: str = "Unknown"
    dfam_order: str = "Unknown"
    dfam_superfamily: str = "Unknown"
    dfam_ct: str = ""


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def nhmmer_tblout_complete(tblout: str | Path) -> bool:
    """Return True only for a tblout that nhmmer finished writing."""
    tblout = Path(tblout)
    if not tblout.exists() or tblout.stat().st_size == 0:
        return False

    last_nonempty = ""
    with open(tblout, errors="replace") as handle:
        for line in handle:
            if line.strip():
                last_nonempty = line.strip()

    return last_nonempty == "# [ok]"


def _tblout_targets_belong_to_chunk(
    tblout: Path,
    chunk_fasta: Path,
) -> bool:
    sequence_lengths: dict[str, int] = {}
    current_id: str | None = None
    with open(chunk_fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                current_id = line[1:].strip().split()[0]
                sequence_lengths[current_id] = 0
            elif current_id is not None:
                sequence_lengths[current_id] += len(line.strip())

    try:
        with open(tblout) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split()
                if len(fields) < 11:
                    return False
                target = fields[0]
                try:
                    target_length = int(fields[10])
                except ValueError:
                    return False
                if sequence_lengths.get(target) != target_length:
                    return False
    except OSError:
        return False

    return True


def _chunk_signature(
    chunk_fasta: Path,
    dfam_db: Path,
    nhmmer_bin: str,
) -> dict:
    return {
        "schema_version": DFAM_CHUNK_CACHE_SCHEMA_VERSION,
        "chunk_sha256": _sha256_file(chunk_fasta),
        "dfam_db": _database_signature(dfam_db),
        "nhmmer_bin": str(nhmmer_bin),
    }


def _database_signature(dfam_db: Path) -> dict:
    db_stat = dfam_db.stat()
    return {
        "path": str(dfam_db.resolve()),
        "size": db_stat.st_size,
        "mtime_ns": db_stat.st_mtime_ns,
    }


def _scan_input_signature(
    consensus_fasta: Path,
    dfam_db: Path,
    nhmmer_bin: str,
) -> dict:
    return {
        "consensus_sha256": _sha256_file(consensus_fasta),
        "dfam_db": _database_signature(dfam_db),
        "nhmmer_bin": str(nhmmer_bin),
    }


def _infer_cached_chunk_count(
    chunks_dir: Path,
    tblout_dir: Path,
) -> int | None:
    tblout_indices = []
    for path in tblout_dir.glob("chunk_*.tblout"):
        match = CHUNK_FILE_PATTERN.match(path.name)
        if match:
            tblout_indices.append(int(match.group(1)))

    if not tblout_indices:
        return None

    checkpoint_indices = list(tblout_indices)
    for path in chunks_dir.glob("chunk_*.fa"):
        match = CHUNK_FILE_PATTERN.match(path.name)
        if match:
            checkpoint_indices.append(int(match.group(1)))

    return max(checkpoint_indices)


def cached_nhmmer_chunk_count(outdir: str | Path) -> int | None:
    outdir = Path(outdir)
    return _infer_cached_chunk_count(
        outdir / "chunks",
        outdir / "tblouts",
    )


def _chunk_manifest_path(tblout: Path) -> Path:
    return Path(f"{tblout}.manifest.json")


def _read_json(path: Path) -> dict | None:
    try:
        with open(path) as handle:
            data = json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None
    return data if isinstance(data, dict) else None


def _write_json_atomic(path: Path, data: dict) -> None:
    temporary = Path(f"{path}.partial")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)

def run_nhmmer(
    dfam_db: str | Path,
    consensus_fasta: str | Path,
    tblout: str | Path,
    nhmmer_bin: str = "nhmmer",
    cpu: int = 1,
) -> Path:
    tblout = Path(tblout)
    tblout.parent.mkdir(parents=True, exist_ok=True)

    log_path = tblout.with_suffix(".log")
    partial_tblout = Path(f"{tblout}.partial")
    partial_log = Path(f"{log_path}.partial")

    cmd = [
        nhmmer_bin,
        "--cpu", str(cpu),
        "--tblout", str(partial_tblout),
        str(dfam_db),
        str(consensus_fasta),
    ]

    try:
        with open(partial_log, "w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)
    except Exception:
        if partial_log.exists():
            os.replace(partial_log, log_path)
        raise

    if not nhmmer_tblout_complete(partial_tblout):
        if partial_log.exists():
            os.replace(partial_log, log_path)
        raise RuntimeError(
            "nhmmer exited without writing a complete tblout: "
            f"{partial_tblout}"
        )

    os.replace(partial_tblout, tblout)
    os.replace(partial_log, log_path)
    return tblout

def parse_nhmmer_tblout(
    tblout: str | Path,
    max_evalue: float = 1e-5,
    min_score: float = 20.0,
    dfam_metadata: dict[str, dict[str, str]] | None = None,
) -> List[DfamHit]:

    hits: List[DfamHit] = []
    dfam_metadata = dfam_metadata or {}

    with open(tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split()

            if len(parts) < 14:
                continue

            family_id = clean_family_id(parts[0])
            model_name = parts[2]
            accession = parts[3]
            ali_start = int(parts[6])
            ali_end = int(parts[7])
            evalue = float(parts[12])
            score = float(parts[13])

            if evalue > max_evalue:
                continue

            if score < min_score:
                continue

            meta = dfam_metadata.get(model_name, {})

            hits.append(
                DfamHit(
                    family_id=family_id,
                    model_name=model_name,
                    accession=accession,
                    evalue=evalue,
                    score=score,
                    ali_start=ali_start,
                    ali_end=ali_end,
                    dfam_class=meta.get("dfam_class", "Unknown"),
                    dfam_order=meta.get("dfam_order", "Unknown"),
                    dfam_superfamily=meta.get(
                        "dfam_superfamily",
                        "Unknown",
                    ),
                    dfam_ct=meta.get("dfam_ct", ""),
                )
            )

    return hits

def best_dfam_hits_by_family(hits: List[DfamHit]) -> dict[str, DfamHit]:
    best = {}

    for hit in hits:
        current = best.get(hit.family_id)

        if current is None:
            best[hit.family_id] = hit
            continue

        if hit.score > current.score:
            best[hit.family_id] = hit
            continue

        if hit.score == current.score and hit.evalue < current.evalue:
            best[hit.family_id] = hit

    return best


def merge_tblout_files(
    tblouts: list[Path],
    merged_tblout: str | Path,
) -> Path:
    merged_tblout = Path(merged_tblout)
    merged_tblout.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(f"{merged_tblout}.partial")

    with open(partial, "w") as out:
        for path in tblouts:
            with open(path) as fh:
                for line in fh:
                    out.write(line)

    os.replace(partial, merged_tblout)
    return merged_tblout

def run_nhmmer_parallel_chunks(
    dfam_db: str | Path,
    consensus_fasta: str | Path,
    outdir: str | Path,
    chunks: int,
    nhmmer_bin: str = "nhmmer",
    cpu_per_job: int = 1,
    max_parallel: int | None = None,
) -> Path:
    """
    Run nhmmer on a consensus FASTA split into length-balanced chunks.

    The merged output is written to:

        <outdir>/dfam.merged.tblout

    Chunk FASTAs are written to:

        <outdir>/chunks/

    Per-chunk tblout files are written to:

        <outdir>/tblouts/
    """

    dfam_db = Path(dfam_db)
    consensus_fasta = Path(consensus_fasta)
    outdir = Path(outdir)

    chunks_dir = outdir / "chunks"
    tblout_dir = outdir / "tblouts"

    chunks_dir.mkdir(parents=True, exist_ok=True)
    tblout_dir.mkdir(parents=True, exist_ok=True)

    requested_chunks = max(1, chunks)
    scan_manifest_path = outdir / "dfam.scan.manifest.json"
    scan_input_signature = _scan_input_signature(
        consensus_fasta,
        dfam_db,
        nhmmer_bin,
    )
    cached_scan_manifest = _read_json(scan_manifest_path)
    scan_manifest_compatible = (
        cached_scan_manifest is not None
        and cached_scan_manifest.get("schema_version")
        == DFAM_SCAN_CACHE_SCHEMA_VERSION
        and cached_scan_manifest.get("input_signature")
        == scan_input_signature
        and isinstance(cached_scan_manifest.get("chunks"), int)
        and cached_scan_manifest["chunks"] >= 1
    )
    inferred_cached_chunks = _infer_cached_chunk_count(
        chunks_dir,
        tblout_dir,
    )
    legacy_checkpoint_detected = (
        not scan_manifest_compatible
        and cached_scan_manifest is None
        and inferred_cached_chunks is not None
    )

    if scan_manifest_compatible:
        effective_chunks = cached_scan_manifest["chunks"]
    elif legacy_checkpoint_detected:
        effective_chunks = inferred_cached_chunks
    else:
        effective_chunks = requested_chunks

    if effective_chunks != requested_chunks:
        print(
            f"[DRayTE][Dfam] cached checkpoint uses "
            f"{effective_chunks} chunks; ignoring requested YAML "
            f"chunks={requested_chunks} while resuming",
            file=sys.stderr,
            flush=True,
        )

    old_chunk_hashes = {
        path.name: _sha256_file(path)
        for path in chunks_dir.glob("chunk_*.fa")
        if path.is_file()
    }

    fasta_chunks = split_fasta_balanced_by_length(
        input_fasta=consensus_fasta,
        outdir=chunks_dir,
        n_chunks=effective_chunks,
        prefix="chunk",
    )

    if not fasta_chunks:
        merged = outdir / "dfam.merged.tblout"
        merged.write_text("")
        return merged

    _write_json_atomic(
        scan_manifest_path,
        {
            "schema_version": DFAM_SCAN_CACHE_SCHEMA_VERSION,
            "chunks": len(fasta_chunks),
            "input_signature": scan_input_signature,
        },
    )

    if max_parallel is None:
        max_parallel = len(fasta_chunks)

    max_parallel = max(1, min(max_parallel, len(fasta_chunks)))

    allow_manifestless_checkpoint = (
        scan_manifest_compatible
        or legacy_checkpoint_detected
    )

    signatures = {
        path: _chunk_signature(path, dfam_db, nhmmer_bin)
        for path in fasta_chunks
    }
    reusable: list[Path] = []
    pending: list[Path] = []

    for chunk_fasta in fasta_chunks:
        tblout = tblout_dir / f"{chunk_fasta.stem}.tblout"
        manifest = _read_json(_chunk_manifest_path(tblout))

        if (
            nhmmer_tblout_complete(tblout)
            and manifest == signatures[chunk_fasta]
        ):
            reusable.append(tblout)
            continue

        old_chunk_hash = old_chunk_hashes.get(chunk_fasta.name)
        legacy_cache_is_compatible = (
            manifest is None
            and allow_manifestless_checkpoint
            and old_chunk_hash is not None
            and nhmmer_tblout_complete(tblout)
            and _tblout_targets_belong_to_chunk(tblout, chunk_fasta)
            and tblout.stat().st_mtime_ns >= dfam_db.stat().st_mtime_ns
        )
        if legacy_cache_is_compatible:
            _write_json_atomic(
                _chunk_manifest_path(tblout),
                signatures[chunk_fasta],
            )
            reusable.append(tblout)
            continue

        pending.append(chunk_fasta)

    print(
        f"[DRayTE][Dfam] resume check: "
        f"{len(reusable)} completed chunks reused, "
        f"{len(pending)} chunks to run",
        file=sys.stderr,
        flush=True,
    )

    def _run_one(chunk_fasta: Path) -> Path:
        tblout = tblout_dir / f"{chunk_fasta.stem}.tblout"

        result = run_nhmmer(
            dfam_db=dfam_db,
            consensus_fasta=chunk_fasta,
            tblout=tblout,
            nhmmer_bin=nhmmer_bin,
            cpu=cpu_per_job,
        )
        _write_json_atomic(
            _chunk_manifest_path(tblout),
            signatures[chunk_fasta],
        )
        return result

    tblouts: list[Path] = list(reusable)

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        futures = {
            executor.submit(_run_one, chunk): chunk
            for chunk in pending
        }

        total = len(fasta_chunks)
        completed = len(reusable)

        for future in as_completed(futures):
            chunk = futures[future]
            tblout = future.result()
            tblouts.append(tblout)

            completed += 1

            print(
                f"[DRayTE][Dfam] progress "
                f"{completed}/{total} chunks completed "
                f"| last={chunk.name}",
                file=sys.stderr,
                flush=True,
            )

    merged = outdir / "dfam.merged.tblout"
    return merge_tblout_files(sorted(tblouts), merged)
