from __future__ import annotations

import fcntl
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List

from Bio import SeqIO

from .homology import homology_from_repeatmasker_header
from .ids import clean_family_id


EXTERNAL_RESCUE_ALGORITHM_VERSION = 2
MMSEQS_OUTPUT_FORMAT = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"


@dataclass
class MMseqsHit:
    query_id: str
    target_id: str
    target_header: str
    identity: float
    aln_len: int
    query_cov: float
    target_cov: float
    qstart: int
    qend: int
    tstart: int
    tend: int
    evalue: float
    bits: float


@dataclass
class RescueEvidence:
    family_id: str
    rescue_class: str
    rescue_order: str
    rescue_superfamily: str
    rescue_identity: float
    rescue_aln_len: int
    rescue_query_cov: float
    rescue_target_cov: float
    rescue_bits: float
    rescue_target: str


def sequence_lengths_by_clean_id(
    fasta: str | Path,
) -> dict[str, int]:
    return {
        clean_family_id(rec.id): len(rec.seq)
        for rec in SeqIO.parse(str(fasta), "fasta")
    }


def run_mmseqs_rescue(
    query_fasta: str | Path,
    reference_fasta: str | Path,
    outdir: str | Path,
    mmseqs_bin: str = "mmseqs",
    threads: int = 1,
    min_seq_id: float = 0.75,
    coverage: float = 0.50,
) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out_tsv = outdir / "mmseqs_rescue.tsv"
    tmpdir = outdir / "tmp"

    cmd = [
        mmseqs_bin,
        "easy-search",
        str(query_fasta),
        str(reference_fasta),
        str(out_tsv),
        str(tmpdir),
        "--search-type",
        "3",
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--cov-mode",
        "0",
        "--threads",
        str(threads),
    ]

    with open(outdir / "mmseqs_rescue.log", "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)

    return out_tsv


def parse_mmseqs_tsv(
    path: str | Path,
    query_lengths: dict[str, int] | None = None,
    target_lengths: dict[str, int] | None = None,
    min_identity: float = 0.80,
    min_aln_len: int = 500,
    min_query_cov: float = 0.0,
    min_target_cov: float = 0.0,
    *,
    skip_self: bool = True,
    strict: bool = False,
) -> List[MMseqsHit]:

    hits: List[MMseqsHit] = []

    query_lengths = query_lengths or {}
    target_lengths = target_lengths or {}

    with open(path) as fh:
        for line_number, line in enumerate(fh, 1):
            if not line.strip():
                continue

            parts = line.rstrip().split("\t")

            if len(parts) < 12:
                if strict:
                    raise ValueError(f"Invalid MMseqs row at {path}:{line_number}: expected 12 columns")
                continue

            query_raw = parts[0]
            target_raw = parts[1]

            query_id = clean_family_id(query_raw)
            target_id = clean_family_id(target_raw)

            if skip_self and query_id == target_id:
                continue

            identity = float(parts[2])
            aln_len = int(parts[3])

            qlen = query_lengths.get(query_id)
            tlen = target_lengths.get(target_id)

            if strict:
                if not qlen or not tlen:
                    raise ValueError(f"MMseqs row contains an unexpected sequence ID at {path}:{line_number}")
                if not math.isfinite(identity) or not 0 <= identity <= 1 or aln_len <= 0:
                    raise ValueError(f"Invalid MMseqs identity/alignment length at {path}:{line_number}")
                if not all(math.isfinite(float(parts[index])) for index in (10, 11)):
                    raise ValueError(f"Non-finite MMseqs score at {path}:{line_number}")

            query_cov = aln_len / qlen if qlen else 0.0
            target_cov = aln_len / tlen if tlen else 0.0

            if identity < min_identity:
                continue

            if aln_len < min_aln_len:
                continue

            if qlen and query_cov < min_query_cov:
                continue

            if tlen and target_cov < min_target_cov:
                continue

            hits.append(
                MMseqsHit(
                    query_id=query_id,
                    target_id=target_id,
                    target_header=target_raw,
                    identity=identity,
                    aln_len=aln_len,
                    query_cov=query_cov,
                    target_cov=target_cov,
                    qstart=int(parts[6]),
                    qend=int(parts[7]),
                    tstart=int(parts[8]),
                    tend=int(parts[9]),
                    evalue=float(parts[10]),
                    bits=float(parts[11]),
                )
            )

    return hits


def best_mmseqs_hits_by_family(
    hits: List[MMseqsHit],
) -> dict[str, MMseqsHit]:
    best = {}

    for hit in hits:
        current = best.get(hit.query_id)

        if current is None:
            best[hit.query_id] = hit
            continue

        better = False

        if hit.bits > current.bits:
            better = True

        elif hit.bits == current.bits:
            if hit.identity > current.identity:
                better = True

            elif (
                hit.identity == current.identity
                and hit.query_cov > current.query_cov
            ):
                better = True

            elif (
                hit.identity == current.identity
                and hit.query_cov == current.query_cov
                and hit.aln_len > current.aln_len
            ):
                better = True

        if better:
            best[hit.query_id] = hit

    return best
def infer_rescue_evidence(
    hits: List[MMseqsHit],
    *,
    min_score_fraction: float = 0.95,
) -> dict[str, RescueEvidence]:
    """Infer consistent labels from qualified one-way similarity hits."""
    evidence, _ = infer_one_way_rescue_evidence(
        hits, min_score_fraction=min_score_fraction,
    )
    return evidence


@dataclass
class ExternalRescueResult:
    forward_hits: list[MMseqsHit]
    evidence: dict[str, RescueEvidence]
    report: dict


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, value: dict):
    temporary = path.with_suffix(path.suffix + "." + uuid.uuid4().hex + ".partial")
    try:
        with temporary.open("x") as handle:
            json.dump(value, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _external_lengths(fasta: Path) -> dict[str, int]:
    lengths = {}
    for record in SeqIO.parse(str(fasta), "fasta"):
        family_id = clean_family_id(record.id)
        if not family_id or family_id in lengths:
            raise ValueError(f"Duplicate or empty clean family ID in MMseqs FASTA {fasta}: {family_id}")
        if not len(record.seq):
            raise ValueError(f"Empty sequence in MMseqs FASTA {fasta}: {record.id}")
        lengths[family_id] = len(record.seq)
    if not lengths:
        raise ValueError(f"MMseqs FASTA contains no records: {fasta}")
    return lengths


def _hit_summary(hit: MMseqsHit) -> dict:
    homology = homology_from_repeatmasker_header(hit.target_header)
    return {
        "target_id": hit.target_id, "target_header": hit.target_header,
        "identity": hit.identity, "aln_len": hit.aln_len,
        "query_cov": hit.query_cov, "target_cov": hit.target_cov,
        "bits": hit.bits, "evalue": hit.evalue,
        "class": homology.homology_class if homology else "Unknown",
        "order": homology.homology_order if homology else "Unknown",
        "superfamily": homology.homology_superfamily if homology else "Unknown",
    }


def infer_one_way_rescue_evidence(
    hits: list[MMseqsHit], *, min_score_fraction: float = 0.95,
) -> tuple[dict[str, RescueEvidence], dict[str, dict]]:
    """Transfer labels supported by competitive query-to-reference matches.

    References may support any number of query families. Alternative reference
    models within ``min_score_fraction`` of the highest bit score must agree on
    class and order. Multiple HSP rows for a reference are never added together.
    """
    if not math.isfinite(min_score_fraction) or not 0 < min_score_fraction <= 1:
        raise ValueError("MMseqs min_score_fraction must be finite and in (0, 1]")
    grouped: dict[str, dict[str, MMseqsHit]] = {}

    def rank(hit):
        # Negative numeric fields put strongest support first; remaining fields
        # make exact ties independent of TSV row order.
        return (-hit.bits, -hit.identity, -hit.query_cov, -hit.aln_len,
                hit.target_id, hit.qstart, hit.qend, hit.tstart, hit.tend,
                hit.evalue)

    for hit in hits:
        if not math.isfinite(hit.bits):
            raise ValueError("MMseqs one-way evidence contains a non-finite bit score")
        if hit.bits <= 0:
            continue
        targets = grouped.setdefault(hit.query_id, {})
        previous = targets.get(hit.target_id)
        if previous is None or rank(hit) < rank(previous):
            targets[hit.target_id] = hit

    evidence = {}
    decisions = {}
    for family_id, by_target in sorted(grouped.items()):
        ordered = sorted(by_target.values(), key=rank)
        best = ordered[0]
        competitive = [hit for hit in ordered
                       if hit.bits >= best.bits * min_score_fraction]
        labels = [homology_from_repeatmasker_header(hit.target_header)
                  for hit in competitive]
        decision = {
            "status": "not_rescued", "selected_target": None,
            "best_target": best.target_header,
            "qualified_targets": len(ordered),
            "near_best_targets": [_hit_summary(hit) for hit in competitive],
        }
        decisions[family_id] = decision
        if any(label is None or label.homology_order == "Unknown" for label in labels):
            decision["reason"] = "unclassified_competing_reference"
            continue
        if len({(label.homology_class, label.homology_order) for label in labels}) != 1:
            decision["reason"] = "conflicting_class_or_order"
            continue
        superfamilies = {label.homology_superfamily for label in labels}
        superfamily = "Unknown"
        if (len(superfamilies) == 1
                and next(iter(superfamilies)) not in {"", "NA", "Unknown", "unknown", None}
                and all(hit.identity >= 0.80 and hit.query_cov >= 0.80
                        and hit.target_cov >= 0.80 for hit in competitive)):
            superfamily = next(iter(superfamilies))
        label = labels[0]
        evidence[family_id] = RescueEvidence(
            family_id=family_id, rescue_class=label.homology_class,
            rescue_order=label.homology_order, rescue_superfamily=superfamily,
            rescue_identity=best.identity, rescue_aln_len=best.aln_len,
            rescue_query_cov=best.query_cov, rescue_target_cov=best.target_cov,
            rescue_bits=best.bits, rescue_target=best.target_header,
        )
        decision.update(status="rescued", selected_target=best.target_header,
                        reason="consistent_reference_support",
                        rescue_class=label.homology_class,
                        rescue_order=label.homology_order,
                        rescue_superfamily=superfamily)
    return evidence, decisions


def _run_reference_mmseqs_rescue(
    query_fasta: str | Path,
    reference_fasta: str | Path,
    outdir: str | Path,
    mmseqs_bin: str = "mmseqs",
    threads: int = 1,
    min_identity: float = 0.80,
    min_aln_len: int = 500,
    min_query_cov: float = 0.40,
    min_target_cov: float = 0.30,
    logger=None,
    *,
    min_score_fraction: float = 0.95,
) -> ExternalRescueResult:
    """Search a classified reference and publish validated query-side rescue.

    Each library has its own ID namespace. The complete checkpoint fingerprints
    both input files, the forward search output and all result-affecting search
    and filtering settings. Unverified legacy TSVs are never adopted.
    """
    if not math.isfinite(min_score_fraction) or not 0 < min_score_fraction <= 1:
        raise ValueError("MMseqs min_score_fraction must be finite and in (0, 1]")
    query = Path(query_fasta).expanduser().resolve()
    reference = Path(reference_fasta).expanduser().resolve()
    directory = Path(outdir).expanduser().absolute()
    if directory.is_symlink():
        raise ValueError(f"Refusing a symlinked MMseqs output directory: {directory}")
    directory.mkdir(parents=True, exist_ok=True)
    forward_path = directory / "mmseqs_rescue.tsv"
    checkpoint = directory / "mmseqs_rescue.checkpoint.json"
    report_path = directory / "mmseqs_rescue.report.json"
    log_path = directory / "mmseqs_rescue.log"
    for path in (forward_path, checkpoint, report_path, log_path):
        if path.is_symlink():
            raise ValueError(f"Refusing a symlinked MMseqs managed output: {path}")
    for name, value in (("min_identity", min_identity), ("min_query_cov", min_query_cov),
                        ("min_target_cov", min_target_cov)):
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"MMseqs {name} must be finite and between 0 and 1")
    if type(min_aln_len) is not int or min_aln_len < 1 or type(threads) is not int or threads < 1:
        raise ValueError("MMseqs min_aln_len and threads must be positive integers")
    query_hash, reference_hash = _sha256(query), _sha256(reference)
    query_lengths = _external_lengths(query)
    reference_lengths = _external_lengths(reference)
    reference_headers = {clean_family_id(record.id): record.id for record in SeqIO.parse(str(reference), "fasta")}
    if _sha256(query) != query_hash or _sha256(reference) != reference_hash:
        raise RuntimeError("MMseqs rescue input changed while reading FASTA metadata")
    search_parameters = {"search_type": 3, "min_seq_id": 0.75, "coverage": 0.50,
                         "cov_mode": 0, "format_output": MMSEQS_OUTPUT_FORMAT}
    filter_parameters = {"min_identity": min_identity, "min_aln_len": min_aln_len,
                         "min_query_cov": min_query_cov, "min_target_cov": min_target_cov,
                         "min_score_fraction": min_score_fraction}
    signature = {"algorithm_version": EXTERNAL_RESCUE_ALGORITHM_VERSION,
                 "mode": "one-way",
                 "query_fasta": str(query), "reference_fasta": str(reference),
                 "query_sha256": query_hash, "reference_sha256": reference_hash,
                 "mmseqs_bin": str(mmseqs_bin), "search_parameters": search_parameters,
                 "filter_parameters": filter_parameters}

    def parse_forward(forward):
        retained = parse_mmseqs_tsv(
            forward, query_lengths, reference_lengths, min_identity=min_identity,
            min_aln_len=min_aln_len, min_query_cov=min_query_cov,
            min_target_cov=min_target_cov, skip_self=False, strict=True,
        )
        raw_forward = parse_mmseqs_tsv(
            forward, query_lengths, reference_lengths, min_identity=0,
            min_aln_len=1, skip_self=False, strict=True,
        )
        # The signed FASTA supplies labels, independent of a tool's header
        # representation or truncation.
        for hit in retained + raw_forward:
            hit.target_header = reference_headers[hit.target_id]
        return retained, raw_forward

    def inputs_unchanged():
        return _sha256(query) == signature["query_sha256"] and _sha256(reference) == signature["reference_sha256"]

    cached = None
    if checkpoint.is_file():
        try:
            value = json.loads(checkpoint.read_text())
            if (value.get("status") == "complete" and value.get("input_signature") == signature
                    and forward_path.is_file()
                    and value.get("forward_sha256") == _sha256(forward_path)):
                cached = parse_forward(forward_path)
        except (OSError, ValueError, TypeError, AttributeError):
            cached = None

    def result_for(hits, reused):
        forward_hits, raw_forward = hits
        evidence, decisions = infer_one_way_rescue_evidence(
            forward_hits, min_score_fraction=min_score_fraction,
        )
        raw_by_query = {}
        for hit in raw_forward:
            raw_by_query.setdefault(hit.query_id, []).append(hit)
        retained_by_query = {}
        for hit in forward_hits:
            retained_by_query.setdefault(hit.query_id, []).append(hit)
        for family_id in sorted(query_lengths):
            raw = raw_by_query.get(family_id, [])
            qualified = retained_by_query.get(family_id, [])
            rejected = []
            for hit in raw:
                failed = []
                for actual, minimum, name in (
                    (hit.identity, min_identity, "identity"),
                    (hit.aln_len, min_aln_len, "alignment_length"),
                    (hit.query_cov, min_query_cov, "query_coverage"),
                    (hit.target_cov, min_target_cov, "reference_coverage"),
                ):
                    if actual < minimum:
                        failed.append(name)
                if hit.bits <= 0:
                    failed.append("nonpositive_bit_score")
                if failed:
                    rejected.append({**_hit_summary(hit), "failed_filters": failed})
            decision = decisions.setdefault(family_id, {
                "status": "not_rescued", "selected_target": None,
                "qualified_targets": len({hit.target_id for hit in qualified}),
                "near_best_targets": [],
                "reason": ("no_search_hit" if not raw else
                           "no_qualified_hit" if not qualified else
                           "no_positive_bit_score"),
            })
            decision.update(raw_hit_rows=len(raw), qualified_hit_rows=len(qualified),
                            rejected_hits=rejected)
        if not set(evidence).issubset(query_lengths):
            raise RuntimeError("External-reference IDs leaked into query rescue evidence")
        report = {**signature, "forward_tsv": str(forward_path),
                  "forward_sha256": _sha256(forward_path),
                  "forward_hits": len(forward_hits), "rescued_families": len(evidence),
                  "families_with_rescue": len(evidence),
                  "cache_reused": reused, "threads": threads, "report_path": str(report_path),
                  "query_families": len(query_lengths), "reference_families": len(reference_lengths),
                  "query_decisions": decisions,
                  "checkpoint": str(checkpoint), "checkpoint_path": str(checkpoint), "log_path": str(log_path)}
        return ExternalRescueResult(forward_hits, evidence, report)

    if cached is not None:
        if not inputs_unchanged():
            raise RuntimeError("MMseqs rescue input changed during cache validation")
        result = result_for(cached, True)
        _atomic_json(report_path, result.report)
        if logger:
            logger.info("Validated one-way MMseqs cache: query=%s reference=%s", query, reference)
        return result

    if logger:
        logger.info("Running one-way MMseqs rescue: query=%s reference=%s", query, reference)
    _atomic_json(checkpoint, {"status": "running", "input_signature": signature})
    published = False
    with tempfile.TemporaryDirectory(prefix=".reference-", dir=directory) as temporary:
        work = Path(temporary)
        backups = {}
        outputs = (forward_path, report_path)
        for index, path in enumerate(outputs):
            backup = work / f"backup-{index}"
            if path.exists():
                shutil.copy2(path, backup)
                backups[path] = backup
        try:
            forward_tmp = work / "forward.tsv"
            with log_path.open("a") as log:
                log.write(f"\nAttempt {datetime.now(timezone.utc).isoformat()}\n")
                command = [str(mmseqs_bin), "easy-search", str(query), str(reference),
                           str(forward_tmp), str(work / "query_to_reference_tmp"),
                           "--search-type", "3", "--min-seq-id", "0.75", "-c", "0.5",
                           "--cov-mode", "0", "--threads", str(threads),
                           "--format-output", MMSEQS_OUTPUT_FORMAT]
                log.write("query_to_reference: " + json.dumps(command) + "\n")
                log.flush()
                subprocess.run(command, check=True, stdout=log, stderr=log)
                if not forward_tmp.is_file():
                    raise RuntimeError(f"MMseqs did not create its query_to_reference result: {forward_tmp}")
            hits = parse_forward(forward_tmp)
            if not inputs_unchanged():
                raise RuntimeError("MMseqs rescue input changed during search")
            # Backups above permit rollback if publication fails after a rename.
            # An incomplete marker can never authorize an unpublished result.
            published = True
            os.replace(forward_tmp, forward_path)
            result = result_for(hits, False)
            _atomic_json(report_path, result.report)
            _atomic_json(checkpoint, {"status": "complete", "input_signature": signature,
                                     "forward_sha256": result.report["forward_sha256"]})
            if logger:
                logger.info("MMseqs rescue (one-way): forward=%d families_with_evidence=%d",
                            len(result.forward_hits), len(result.evidence))
            return result
        except BaseException as error:
            if published:
                for path in outputs:
                    if path in backups:
                        os.replace(backups[path], path)
                    else:
                        path.unlink(missing_ok=True)
            try:
                _atomic_json(checkpoint, {"status": "failed", "input_signature": signature,
                                         "error": f"{type(error).__name__}: {error}"})
            except OSError:
                pass  # The previously persisted running marker remains incomplete.
            raise


def run_reference_mmseqs_rescue(
    query_fasta: str | Path,
    reference_fasta: str | Path,
    outdir: str | Path,
    mmseqs_bin: str = "mmseqs",
    threads: int = 1,
    min_identity: float = 0.80,
    min_aln_len: int = 500,
    min_query_cov: float = 0.40,
    min_target_cov: float = 0.30,
    logger=None,
    *,
    min_score_fraction: float = 0.95,
) -> ExternalRescueResult:
    """Run or reuse signed one-way reference rescue under a lock."""
    directory = Path(outdir).expanduser().absolute()
    if directory.is_symlink():
        raise ValueError(f"Refusing a symlinked MMseqs output directory: {directory}")
    directory.mkdir(parents=True, exist_ok=True)
    lock = directory / "mmseqs_rescue.lock"
    if lock.is_symlink():
        raise ValueError(f"Refusing a symlinked MMseqs lock: {lock}")
    with lock.open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"Another MMseqs rescue process holds the output lock: {directory}") from error
        try:
            return _run_reference_mmseqs_rescue(
                query_fasta, reference_fasta, directory, mmseqs_bin, threads,
                min_identity, min_aln_len, min_query_cov, min_target_cov, logger,
                min_score_fraction=min_score_fraction,
            )
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)
