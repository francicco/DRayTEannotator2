"""One-time recovery of a legacy classified library from completed curation.

Preview is read-only. Recovery reconstructs scientific FASTA records, rather
than claiming byte-identical restoration or a new RepeatClassifier execution.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import shutil
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from drayte.curation.finalize import normalize_record_header_for_repeatmasker
from drayte.pipeline.config import load_config
from drayte.pipeline.reclassify import (
    checkpoint_data, checkpoint_is_valid, classified_output_has_complete_family_set, write_json_atomic,
)
from drayte.storage.policy import output_lock


class RecoveryError(ValueError):
    pass


def _digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _state(path):
    info = path.stat()
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


def _real_path(root, path, *, directory=False, allow_empty=False):
    path = Path(path)
    if not path.is_absolute():
        raise RecoveryError(f"Recovery proof must use an absolute path: {path}")
    try:
        relative = path.relative_to(root)
    except ValueError as error:
        raise RecoveryError(f"Recovery input is outside the run: {path}") from error
    if ".." in relative.parts:
        raise RecoveryError(f"Noncanonical recovery path: {path}")
    current = root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise RecoveryError(f"Recovery refuses symlinked input: {current}")
    if not (path.is_dir() if directory else path.is_file()):
        raise RecoveryError(f"Required recovery {'directory' if directory else 'file'} is missing: {path}")
    if not directory and not allow_empty and not path.stat().st_size:
        raise RecoveryError(f"Empty recovery input: {path}")
    return path


@dataclass
class Inspection:
    report: dict
    snapshots: dict
    absent: set
    records: list


def _records(path):
    values = {}
    for record in SeqIO.parse(str(path), "fasta"):
        family = record.id.split("#", 1)[0]
        if not family or family in values or not len(record.seq):
            raise RecoveryError(f"Duplicate/empty FASTA family in {path}: {family}")
        values[family] = record
    if not values:
        raise RecoveryError(f"No FASTA records in {path}")
    return values


def _writer_locked(root):
    lock = root / ".storage/run.lock"
    if not lock.exists():
        return False
    _real_path(root, lock, allow_empty=True)
    with lock.open("r") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        finally:
            try:
                fcntl.flock(handle, fcntl.LOCK_UN)
            except OSError:
                pass
    return False


def _inspect(config, *, check_lock=True):
    requested = config.outdir_path.expanduser().absolute()
    if requested.is_symlink() or not requested.is_dir():
        raise RecoveryError(f"Output root must be a real existing directory: {requested}")
    root = requested.resolve()
    species = config.species
    merged = root / f"reclassify/{species}_extended_rep.fa"
    output = Path(str(merged) + ".classified")
    checkpoint = root / "reclassify/repeatclassifier.checkpoint.json"
    receipt = root / "reclassify/reclassify.manifest.json"
    curated = root / "curation/Final.RepeatModeler.Lib.fa"
    curation_receipt = root / "curation/curation.manifest.json"
    sources_dir = root / "extension/final_consensuses"
    snapshots = {}
    absent = set()

    def capture(path, *, directory=False, allow_empty=False):
        path = _real_path(root, path, directory=directory, allow_empty=allow_empty)
        snapshots[path] = _state(path)
        if directory:
            return path, {}
        data = {"path": str(path), "sha256": _digest(path),
                "bytes": snapshots[path][2], "mtime_ns": snapshots[path][3],
                "ctime_ns": snapshots[path][4]}
        if _state(path) != snapshots[path]:
            raise RecoveryError(f"Recovery input changed while being read: {path}")
        return path, data

    def load_receipt(path):
        path, info = capture(path)
        value = json.loads(path.read_text())
        if not isinstance(value, dict):
            raise RecoveryError(f"Recovery receipt is not a JSON mapping: {path}")
        return value, info

    def points(value, expected, label):
        if not isinstance(value, str) or not Path(value).is_absolute() or Path(value).resolve() != expected.resolve():
            raise RecoveryError(f"{label} does not identify the expected file: {expected}")

    merged, input_info = capture(merged)
    raw = _records(merged)
    binary = str(config.extra.get("repeatclassifier_bin", "RepeatClassifier"))
    report = {"read_only": True, "eligible": False, "errors": [], "outdir": str(root),
              "current_input": input_info, "family_count": len(raw), "applied": False}
    if output.exists() or output.is_symlink():
        _, report["current_output"] = capture(output, allow_empty=True)
    else:
        absent.add(output)
    if checkpoint.exists() or checkpoint.is_symlink():
        marker, marker_info = load_receipt(checkpoint)
        report["current_checkpoint"] = {**marker_info, "status": marker.get("status")}
        if marker.get("status") == "complete":
            if (checkpoint_is_valid(checkpoint, merged, output, binary)
                    and classified_output_has_complete_family_set(merged, output)):
                report.update(already_complete=True, explanation="The modern completion checkpoint is already valid")
                return Inspection(report, snapshots, absent, [])
            raise RecoveryError("A modern completed checkpoint exists but does not match; legacy recovery cannot bypass it")
        if marker.get("status") != "running" or marker.get("schema_version") != 1:
            raise RecoveryError("Only an absent checkpoint or a matching abandoned running attempt can be recovered")
        points(marker.get("input_fasta"), merged, "Running checkpoint input")
        points(marker.get("output_fasta"), output, "Running checkpoint output")
        if marker.get("input_sha256") != input_info["sha256"] or marker.get("repeatclassifier_bin") != binary:
            raise RecoveryError("Running checkpoint input hash or executable does not match the requested recovery")
    else:
        absent.add(checkpoint)

    old, old_info = load_receipt(receipt)
    if old.get("stage") != "reclassify":
        raise RecoveryError("The historical RepeatClassifier manifest does not record a completed reclassify stage")
    points(old.get("outdir"), root / "reclassify", "Historical reclassify directory")
    points(old.get("merged_library"), merged, "Historical reclassify input")
    points(old.get("classified_library"), output, "Historical reclassify output")
    if "input_sha256" in old and old["input_sha256"] != input_info["sha256"]:
        raise RecoveryError("The historical reclassify receipt records a different input hash")

    curated_manifest, curation_info = load_receipt(curation_receipt)
    if curated_manifest.get("stage") != "curation":
        raise RecoveryError("The curation manifest does not record a completed curation stage")
    points(curated_manifest.get("outdir"), root / "curation", "Completed curation directory")
    library_key = "final_library" if "final_library" in curated_manifest else "curated_library"
    points(curated_manifest.get(library_key), curated, "Completed curation library")
    if "curated_library" in curated_manifest:
        points(curated_manifest["curated_library"], curated, "Completed curated library")
    curated, curated_info = capture(curated)
    table_path = Path(curated_manifest.get("family_table", ""))
    table_path, table_info = capture(table_path)
    if old_info["mtime_ns"] > curation_info["mtime_ns"] or old_info["ctime_ns"] > curation_info["mtime_ns"]:
        raise RecoveryError("The historical reclassification receipt does not predate completed curation")
    for label, info in (("curated library", curated_info), ("family table", table_info)):
        if max(info["mtime_ns"], info["ctime_ns"]) > curation_info["mtime_ns"]:
            raise RecoveryError(f"The {label} changed after its curation completion receipt")

    capture(sources_dir, directory=True)
    sources = sorted(sources_dir.glob(f"{species}*.fa")) or sorted(sources_dir.glob("*.fa"))
    if not sources:
        raise RecoveryError("No extension consensus sources were found")
    if "n_input_fastas" in old and old["n_input_fastas"] != len(sources):
        raise RecoveryError("Historical input FASTA count does not match the extension sources")
    concatenated = hashlib.sha256()
    newest_mtime = newest_ctime = 0
    for source in sources:
        source, info = capture(source)
        newest_mtime = max(newest_mtime, info["mtime_ns"])
        newest_ctime = max(newest_ctime, info["ctime_ns"])
        if max(info["mtime_ns"], info["ctime_ns"]) > old_info["mtime_ns"]:
            raise RecoveryError(f"An extension source changed after historical RepeatClassifier completion: {source}")
        with source.open(encoding="utf-8") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), ""):
                concatenated.update(block.encode("utf-8"))
    if concatenated.hexdigest() != input_info["sha256"]:
        raise RecoveryError("Current merged input does not match the extension consensus sources")

    candidate = _records(curated)
    if candidate.keys() != raw.keys():
        raise RecoveryError("The curated library does not preserve the exact input family set")
    counts = {"identical_sequences": 0, "reverse_complement_sequences": 0}
    for family, record in candidate.items():
        first, second = str(raw[family].seq).upper(), str(record.seq).upper()
        if len(first) != len(second):
            raise RecoveryError(f"Curated family length differs from current input: {family}")
        if first == second:
            counts["identical_sequences"] += 1
        elif first == str(record.seq.reverse_complement()).upper():
            counts["reverse_complement_sequences"] += 1
        else:
            raise RecoveryError(f"Curated family contains differences other than orientation: {family}")

    headers = {}
    with table_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"original_header", "name", "consensus_len"}.issubset(reader.fieldnames or []):
            raise RecoveryError("The family table lacks original headers, names, or consensus lengths")
        for row in reader:
            family, header = row["name"], row["original_header"]
            if family in headers or family not in candidate or header.split("#", 1)[0] != family:
                raise RecoveryError(f"Invalid or duplicate family-table record: {family}")
            if int(row["consensus_len"]) != len(candidate[family].seq):
                raise RecoveryError(f"Family-table length differs from curated sequence: {family}")
            reconstructed = SeqRecord(candidate[family].seq, id=header, description=header)
            if normalize_record_header_for_repeatmasker(reconstructed).id != candidate[family].id:
                raise RecoveryError(f"Curated classification is inconsistent with its historical header: {family}")
            headers[family] = header
    if headers.keys() != candidate.keys():
        raise RecoveryError("The family table does not cover the exact curated family set")
    recovered = [SeqRecord(record.seq, id=headers[family], description=headers[family])
                 for family, record in candidate.items()]
    report.update(eligible=True, recovered_records=len(recovered), sequence_comparison=counts,
                  historical_reclassify_receipt=old_info, curation_receipt=curation_info,
                  curated_source=curated_info, historical_family_table=table_info,
                  extension_source_proof={"files": len(sources), "newest_mtime_ns": newest_mtime,
                                          "newest_ctime_ns": newest_ctime,
                                          "concatenated_sha256": concatenated.hexdigest()},
                  reconstruction="Original classification IDs from the historical family table; oriented sequence bodies from completed curation. Not byte-identical FASTA restoration or a fresh classifier run.")
    _assert_unchanged(snapshots, absent)
    if check_lock and _writer_locked(root):
        report.update(eligible=False, errors=["An active updated pipeline/storage process holds the output lock; stop or finish it before recovery"])
    return Inspection(report, snapshots, absent, recovered)


def _assert_unchanged(snapshots, absent):
    for path, previous in snapshots.items():
        if path.is_symlink() or not path.exists() or _state(path) != previous:
            raise RecoveryError(f"Recovery input changed during inspection: {path}")
    for path in absent:
        if path.exists() or path.is_symlink():
            raise RecoveryError(f"A previously absent recovery output appeared during inspection: {path}")


def audit_recovery(config):
    """Inspect receipts and source equivalence without creating any files."""
    try:
        return _inspect(config).report
    except (OSError, ValueError, KeyError, TypeError) as error:
        return {"read_only": True, "eligible": False, "applied": False,
                "outdir": str(config.outdir_path), "errors": [str(error)]}


def apply_recovery(config):
    """Reconstruct only after revalidation under the pipeline's output lock."""
    root = config.outdir_path.expanduser().resolve()
    with output_lock(root):
        inspected = _inspect(config, check_lock=False)
        if inspected.report.get("already_complete"):
            return {**inspected.report, "read_only": False}
        if not inspected.report["eligible"]:
            raise RecoveryError("Recovery preconditions did not pass")
        directory = root / "reclassify"
        merged = directory / f"{config.species}_extended_rep.fa"
        output = Path(str(merged) + ".classified")
        checkpoint = directory / "repeatclassifier.checkpoint.json"
        receipt = directory / "reclassify.manifest.json"
        sidecar = directory / "repeatclassifier.recovery.json"
        if sidecar.exists() or sidecar.is_symlink():
            _real_path(root, sidecar, allow_empty=True)
        backups = directory / "recovery_backups"
        if backups.is_symlink():
            raise RecoveryError("Recovery backup directory cannot be a symlink")
        identity = uuid.uuid4().hex
        temporary = directory / f".recovered-{identity}.fa.partial"
        mutable = (output, checkpoint, receipt, sidecar)
        backup_dir = backups / identity
        saved = {}
        initial_mutable_state = {path: _state(path) if path.exists() else None for path in mutable}
        try:
            with temporary.open("x") as handle:
                SeqIO.write(inspected.records, handle, "fasta")
                handle.flush()
                os.fsync(handle.fileno())
            written = _records(temporary)
            if written.keys() != _records(merged).keys():
                raise RecoveryError("Reconstructed FASTA family set changed unexpectedly")
            for intended in inspected.records:
                actual = written[intended.id.split("#", 1)[0]]
                if actual.id != intended.id or str(actual.seq) != str(intended.seq):
                    raise RecoveryError("Reconstructed FASTA content differs from its validated source")
            _assert_unchanged(inspected.snapshots, inspected.absent)
            backup_dir.mkdir(parents=True, exist_ok=False)
            for path in mutable:
                if path.exists():
                    backup = backup_dir / path.name
                    shutil.copy2(path, backup)
                    if _digest(path) != _digest(backup):
                        raise RecoveryError(f"Backup validation failed: {path}")
                    saved[path] = backup
            _assert_unchanged(inspected.snapshots, inspected.absent)
            report = {**inspected.report, "read_only": False, "applied": True,
                      "recovered_at_utc": datetime.now(timezone.utc).isoformat(),
                      "backup_directory": str(backup_dir), "fresh_classifier_execution": False}
            mutated = False
            try:
                os.replace(temporary, output)
                mutated = True
                data = checkpoint_data(merged, output, str(config.extra.get("repeatclassifier_bin", "RepeatClassifier")))
                report["recovered_output"] = {"path": str(output), "sha256": data["output_sha256"]}
                # Publish reconstruction provenance before a normal complete
                # checkpoint can become visible, including across a hard kill.
                write_json_atomic(sidecar, report)
                write_json_atomic(checkpoint, data)
                write_json_atomic(receipt, {"stage": "reclassify", "outdir": str(directory),
                                           "merged_library": str(merged), "classified_library": str(output),
                                           "n_input_fastas": inspected.report["extension_source_proof"]["files"],
                                           "input_sha256": data["input_sha256"], "output_sha256": data["output_sha256"],
                                           "recovery_provenance": str(sidecar)})
            except BaseException:
                if mutated:
                    # Existing allocations are moved back, so rollback does not
                    # need a second full FASTA allocation on a full filesystem.
                    for path in mutable:
                        if path.exists() and _state(path) == initial_mutable_state[path]:
                            continue
                        if path in saved:
                            os.replace(saved[path], path)
                        elif path.exists():
                            path.unlink()
                raise
            return report
        finally:
            temporary.unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--apply", action="store_true", help="Recover an idle run after all provenance checks pass")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    try:
        report = apply_recovery(config) if args.apply else audit_recovery(config)
        print(json.dumps(report, indent=2))
        return 0 if report.get("eligible") or report.get("already_complete") else 1
    except (OSError, ValueError, RuntimeError, KeyError, TypeError) as error:
        print(json.dumps({"read_only": not args.apply, "eligible": False, "applied": False,
                          "errors": [str(error)]}, indent=2))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
