"""Content-based completion records for RepeatMasker stages.

A record is written as ``running`` before execution and as ``complete`` only
after the successful process has produced every required output. A stale table
alone can therefore never turn a failed managed attempt into a resumed stage.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


SCHEMA_VERSION = 1


def checkpoint_path(outdir: Path) -> Path:
    return Path(outdir) / "repeatmasker.checkpoint.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_signature(genome: Path, library: Path, options: list[str]) -> dict:
    return {
        "genome_sha256": sha256_file(genome),
        "library_sha256": sha256_file(library),
        "options": [str(value) for value in options],
    }


def _write(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def begin_attempt(path: Path, signature: dict) -> None:
    _write(path, {"schema_version": SCHEMA_VERSION, "state": "running", "input_signature": signature})


def validate_outputs(outputs: list[Path]) -> None:
    missing = [str(path) for path in outputs if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError("RepeatMasker did not produce required nonempty outputs: " + ", ".join(missing))


def write_checkpoint(path: Path, signature: dict, outputs: list[Path], *, legacy=False) -> None:
    validate_outputs(outputs)
    _write(path, {
        "schema_version": SCHEMA_VERSION,
        "state": "complete",
        "input_signature": signature,
        "outputs": {output.name: {"sha256": sha256_file(output), "size": output.stat().st_size}
                    for output in outputs},
        "legacy_adopted": bool(legacy),
    })


def checkpoint_is_valid(path: Path, signature: dict, outputs: list[Path]) -> bool:
    try:
        record = json.loads(path.read_text())
        if (record.get("schema_version") != SCHEMA_VERSION or record.get("state") != "complete"
                or record.get("input_signature") != signature):
            return False
        validate_outputs(outputs)
        expected = {output.name: {"sha256": sha256_file(output), "size": output.stat().st_size}
                    for output in outputs}
        return record.get("outputs") == expected
    except (OSError, ValueError, RuntimeError, AttributeError):
        return False


def legacy_outputs_usable(outputs: list[Path], log_file: Path) -> bool:
    """Accept complete legacy outputs unless a logged latest attempt failed.

Legacy records have no historical input hashes. This is compatibility adoption,
not proof of historical tool/database provenance; adopted records say so.
"""
    try:
        validate_outputs(outputs)
        latest = None
        if log_file.exists():
            with log_file.open(errors="replace") as handle:
                for line in handle:
                    if line.startswith("=== RepeatMasker started "):
                        latest = "running"
                    elif line.startswith("=== RepeatMasker completed successfully ==="):
                        latest = "complete"
                    elif line.startswith("=== RepeatMasker failed "):
                        latest = "failed"
        return latest in (None, "complete")
    except (OSError, RuntimeError):
        return False


def remove_stale_source_outputs(genome: Path, suffixes: list[str]) -> None:
    """Remove only previous attempt products beside the input genome.

Published stage outputs remain in their stage directory until successful new
outputs are moved into place. The input FASTA is never removed.
"""
    for suffix in suffixes:
        path = Path(str(genome) + suffix)
        if path.exists() or path.is_symlink():
            if path.is_dir() and not path.is_symlink():
                raise RuntimeError(f"Refusing directory where RepeatMasker output was expected: {path}")
            path.unlink()
