#!/usr/bin/env python3
"""Standalone, reversible compaction of completed EarlGrey runs.

Copy this single file anywhere and run it with Python 3.10+ on Linux or macOS.
Only the Python standard library is required; no DrayAnn installation is needed.

Preview: python3 compact_earlgrey.py --root Avcr_Earlgrey
Apply:   python3 compact_earlgrey.py --root Avcr_Earlgrey --apply --idle
Restore: python3 compact_earlgrey.py --root Avcr_Earlgrey --restore --apply --idle

Native EarlGrey does not understand managed archives. Restore before invoking
it directly, or use --run to restore, execute EarlGrey, then compact on success.
The --idle flag asserts no other jobs or workers use the selected output tree.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import logging
import os
import re
import shutil
import signal
import stat
import subprocess
import tarfile
import threading
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

# Self-contained archive engine. Schema 1 remains compatible with archives
# created by the earlier packaged utility; checksums and restore guards apply
# equally to both. This file must remain runnable when copied on its own.

def _sha256(handle) -> str:
    digest = hashlib.sha256()
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
    return digest.hexdigest()


def _file_sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return _sha256(handle)


def _present(path: Path) -> bool:
    return path.exists() or path.is_symlink()


class ArchiveStore:
    """Archive files/trees below ``root`` and restore their original paths."""

    def __init__(self, root: Path):
        self.root = Path(root).expanduser().absolute()
        if self.root.is_symlink() or not self.root.is_dir():
            raise ValueError(f"Archive root must be an existing real directory: {self.root}")
        self.storage = self.root / ".storage"
        self.directory = self.storage / "archives"
        self.manifest = self.storage / "archives.json"
        self._check_ancestors(self.manifest)

    @staticmethod
    def _relative(value, *, storage=False) -> str:
        text = str(value)
        path = PurePosixPath(text)
        if (not text or path.is_absolute() or ".." in path.parts or "\\" in text
                or not path.parts or path == PurePosixPath(".")):
            raise ValueError(f"Unsafe archive path: {value!r}")
        if not storage and path.parts[0] == ".storage":
            raise ValueError("The storage catalog cannot archive itself")
        return path.as_posix()

    def _check_ancestors(self, path: Path):
        try:
            relative = path.relative_to(self.root)
        except ValueError as error:
            raise ValueError(f"Path outside archive root: {path}") from error
        current = self.root
        for part in relative.parts[:-1]:
            current = current / part
            if current.is_symlink() or (current.exists() and not current.is_dir()):
                raise ValueError(f"Unsafe archive ancestor: {current}")

    def _path(self, relative: str) -> Path:
        path = self.root / self._relative(relative)
        self._check_ancestors(path)
        return path

    def _check_link(self, name: str, target: str):
        if not isinstance(target, str) or not target:
            raise ValueError(f"Invalid symlink target in archive: {name}")
        link = self.root / name
        try:
            destination = (link.parent / target).resolve()
            destination.relative_to(self.root.resolve())
        except (ValueError, RuntimeError) as error:
            raise ValueError(f"Symlink escapes archive root or loops: {name} -> {target}") from error

    def _inventory(self, relative: str) -> list[dict]:
        source = self._path(relative)
        if source.is_symlink():
            raise ValueError(f"Cannot archive a symlink as the source root: {source}")
        members = []

        def visit(path):
            info = path.lstat()
            name = path.relative_to(self.root).as_posix()
            member = {"name": name, "mode": stat.S_IMODE(info.st_mode), "mtime_ns": info.st_mtime_ns}
            if stat.S_ISDIR(info.st_mode):
                member["type"] = "directory"
                members.append(member)
                for child in sorted(path.iterdir()):
                    visit(child)
            elif stat.S_ISREG(info.st_mode):
                member.update(type="file", size=info.st_size, sha256=_file_sha256(path))
                after = path.lstat()
                if (info.st_ino, info.st_size, info.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
                    raise RuntimeError(f"Source changed while hashing: {path}")
                members.append(member)
            elif stat.S_ISLNK(info.st_mode):
                target = os.readlink(path)
                self._check_link(name, target)
                member.update(type="symlink", linkname=target)
                members.append(member)
            else:
                raise ValueError(f"Unsupported special file in archive source: {path}")

        visit(source)
        return members

    def _validate_record(self, relative: str, record: dict):
        if not isinstance(record, dict):
            raise ValueError(f"Invalid archive record for {relative}")
        archive_name = self._relative(record.get("archive", ""), storage=True)
        if not archive_name.startswith(".storage/archives/") or len(PurePosixPath(archive_name).parts) != 3:
            raise ValueError(f"Archive payload is outside its storage directory: {archive_name}")
        path = self.root / archive_name
        self._check_ancestors(path)
        if path.is_symlink() or not path.is_file():
            raise FileNotFoundError(f"Managed archive is missing or symlinked: {path}")
        if not isinstance(record.get("members"), list) or not record["members"]:
            raise ValueError(f"Archive has no member inventory: {relative}")
        if not self._is_digest(record.get("archive_sha256")):
            raise ValueError(f"Archive has no valid payload checksum: {relative}")
        names = set()
        for member in record["members"]:
            if not isinstance(member, dict) or not isinstance(member.get("name"), str):
                raise ValueError(f"Invalid archive member inventory: {relative}")
            name = self._relative(member["name"])
            if name != member["name"]:
                raise ValueError(f"Noncanonical archive member name: {member['name']}")
            if name != relative and not name.startswith(relative + "/"):
                raise ValueError(f"Archive member is outside its source: {name}")
            if name in names or member.get("type") not in {"file", "directory", "symlink"}:
                raise ValueError(f"Duplicate or invalid archive member: {name}")
            names.add(name)
            if (not isinstance(member.get("mode"), int) or not 0 <= member["mode"] <= 0o7777
                    or not isinstance(member.get("mtime_ns"), int)):
                raise ValueError(f"Invalid file metadata in archive: {name}")
            if member["type"] == "file" and (
                not isinstance(member.get("size"), int) or member["size"] < 0
                or not self._is_digest(member.get("sha256"))
            ):
                raise ValueError(f"Invalid file checksum/size in archive: {name}")
            if member["type"] == "symlink":
                self._check_link(name, member["linkname"])
        if relative not in names:
            raise ValueError(f"Archive inventory is missing its source root: {relative}")
        types = {member["name"]: member["type"] for member in record["members"]}
        if types[relative] == "symlink":
            raise ValueError("Managed source roots cannot be symlinks")
        for name in names:
            for parent in PurePosixPath(name).parents:
                if parent.as_posix() in types and types[parent.as_posix()] != "directory":
                    raise ValueError(f"Archive member has a non-directory parent: {name}")
        return path

    @staticmethod
    def _is_digest(value):
        return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)

    def records(self) -> dict:
        """Return catalog records, failing explicitly on missing payloads."""
        self._check_ancestors(self.manifest)
        if self.manifest.is_symlink():
            raise ValueError("Archive manifest cannot be a symlink")
        if not self.manifest.exists():
            return {}
        data = json.loads(self.manifest.read_text())
        if not isinstance(data, dict) or data.get("schema_version") != 1 or not isinstance(data.get("archives"), dict):
            raise ValueError("Unsupported or invalid archive manifest")
        records = data["archives"]
        for relative, record in records.items():
            if self._relative(relative) != relative:
                raise ValueError(f"Noncanonical archive source: {relative}")
            self._validate_record(relative, record)
        for name in records:
            if any(parent.as_posix() in records for parent in PurePosixPath(name).parents):
                raise ValueError("Nested managed archives are not supported")
        return records

    def _write_records(self, records):
        self._check_ancestors(self.manifest)
        self.directory.mkdir(parents=True, exist_ok=True)
        temporary = self.storage / f"archives.json.{uuid.uuid4().hex}.partial"
        try:
            with temporary.open("x") as handle:
                json.dump({"schema_version": 1, "archives": records}, handle, indent=2, sort_keys=True)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, self.manifest)
            self._sync_directory(self.storage)
        finally:
            temporary.unlink(missing_ok=True)

    def _verify_archive(self, record):
        path = self.root / record["archive"]
        if _file_sha256(path) != record["archive_sha256"]:
            raise ValueError(f"Archive checksum mismatch: {path}")
        expected = {member["name"]: member for member in record["members"]}
        seen = set()
        with tarfile.open(path, "r:gz") as archive:
            for member in archive:
                name = self._relative(member.name)
                if name != member.name or name not in expected or name in seen:
                    raise ValueError(f"Unexpected or duplicate archive member: {name}")
                seen.add(name)
                info = expected[name]
                if info["type"] == "file":
                    if not member.isfile() or member.size != info["size"]:
                        raise ValueError(f"Archive member type/size mismatch: {name}")
                    with archive.extractfile(member) as handle:
                        if _sha256(handle) != info["sha256"]:
                            raise ValueError(f"Archive member checksum mismatch: {name}")
                elif info["type"] == "directory":
                    if not member.isdir():
                        raise ValueError(f"Archive directory type mismatch: {name}")
                elif not member.issym() or member.linkname != info["linkname"]:
                    raise ValueError(f"Archive symlink mismatch: {name}")
        if seen != set(expected):
            raise ValueError(f"Archive is missing members: {path}")

    @staticmethod
    def _same_content(first, second):
        keys = {"name", "type", "size", "sha256", "linkname"}
        return {key: first.get(key) for key in keys} == {key: second.get(key) for key in keys}

    @staticmethod
    def _sync_directory(path):
        descriptor = os.open(path, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)

    def archive(self, relative: str | Path, *, retain=()) -> dict:
        """Archive a complete tree, optionally leaving selected regular files loose.

        Retained names are relative to the store root. The full inventory stays
        in the archive, so restoration also checks the retained files before
        filling in missing scratch. This permits downstream readers to work
        without expanding a completed family's temporary files.
        """
        relative = self._relative(relative)
        retain = {self._relative(name) for name in retain}
        source = self._path(relative)
        records = self.records()
        if relative in records:
            self._verify_archive(records[relative])
            if not _present(source):
                return records[relative]
            self._restore_one(relative, records)
            records = self.records()
        for name in records:
            if name.startswith(relative + "/") or relative.startswith(name + "/"):
                raise ValueError(f"Source overlaps an existing archive: {relative}")
        members = self._inventory(relative)
        by_name = {member["name"]: member for member in members}
        for name in retain:
            if name not in by_name or by_name[name]["type"] != "file":
                raise ValueError(f"Retained archive member must be a regular file: {name}")
        self.directory.mkdir(parents=True, exist_ok=True)
        self._check_ancestors(self.directory / "payload")
        filename = uuid.uuid4().hex + ".tar.gz"
        payload = self.directory / filename
        temporary = self.directory / (filename + ".partial")
        record = {"archive": payload.relative_to(self.root).as_posix(), "members": members}
        try:
            with tarfile.open(temporary, "w:gz", compresslevel=6) as archive:
                for info in members:
                    member = tarfile.TarInfo(info["name"])
                    member.mode = info["mode"]
                    member.mtime = info["mtime_ns"] / 1_000_000_000
                    if info["type"] == "directory":
                        member.type = tarfile.DIRTYPE
                        archive.addfile(member)
                    elif info["type"] == "symlink":
                        member.type = tarfile.SYMTYPE
                        member.linkname = info["linkname"]
                        archive.addfile(member)
                    else:
                        member.size = info["size"]
                        with (self.root / info["name"]).open("rb") as handle:
                            archive.addfile(member, handle)
            with temporary.open("rb") as handle:
                os.fsync(handle.fileno())
            os.replace(temporary, payload)
            self._sync_directory(self.directory)
            record["archive_sha256"] = _file_sha256(payload)
            self._verify_archive(record)
            if self._inventory(relative) != members:
                raise RuntimeError(f"Source changed while creating archive: {source}")
            records[relative] = record
            self._write_records(records)
        except BaseException:
            temporary.unlink(missing_ok=True)
            # Before commit, source remains authoritative. Never delete a
            # payload after its recovery record has been committed.
            try:
                committed = self.records().get(relative, {}).get("archive") == record["archive"]
            except Exception:
                committed = True
            if not committed:
                payload.unlink(missing_ok=True)
            raise
        if retain:
            keep = set(retain)
            for name in retain:
                keep.update(parent.as_posix() for parent in PurePosixPath(name).parents)
            # Delete from the verified inventory, deepest entries first. Never
            # follow symlinks or recursively remove a retained parent directory.
            for info in reversed(members):
                if info["name"] in keep:
                    continue
                path = self._path(info["name"])
                if info["type"] == "directory":
                    path.rmdir()
                else:
                    path.unlink()
        elif source.is_dir():
            if not shutil.rmtree.avoids_symlink_attacks:
                raise RuntimeError("Python lacks symlink-resistant directory removal")
            shutil.rmtree(source)
        else:
            source.unlink()
        return record

    def _restore_one(self, relative, records):
        record = records[relative]
        self._validate_record(relative, record)
        self._verify_archive(record)
        expected = {member["name"]: member for member in record["members"]}
        source = self._path(relative)
        if _present(source):
            for member in self._inventory(relative):
                if member["name"] not in expected or not self._same_content(member, expected[member["name"]]):
                    raise ValueError(f"Conflicting existing file during restore: {member['name']}")
        for info in record["members"]:
            self._check_ancestors(self.root / info["name"])
        staging = self.storage / ("restore-" + uuid.uuid4().hex)
        staging.mkdir()
        try:
            with tarfile.open(self.root / record["archive"], "r:gz") as archive:
                for member in archive:
                    info = expected[member.name]
                    destination = self.root / member.name
                    self._check_ancestors(destination)
                    if info["type"] == "directory":
                        destination.mkdir(parents=True, exist_ok=True)
                    elif not _present(destination):
                        destination.parent.mkdir(parents=True, exist_ok=True)
                        temporary = staging / uuid.uuid4().hex
                        if info["type"] == "file":
                            with archive.extractfile(member) as inp, temporary.open("xb") as out:
                                shutil.copyfileobj(inp, out, 1024 * 1024)
                        else:
                            self._check_link(member.name, info["linkname"])
                            temporary.symlink_to(info["linkname"])
                        os.replace(temporary, destination)
            for info in reversed(record["members"]):
                path = self.root / info["name"]
                if info["type"] != "symlink":
                    os.chmod(path, info["mode"])
                os.utime(path, ns=(info["mtime_ns"], info["mtime_ns"]), follow_symlinks=False)
            if self._inventory(relative) != record["members"]:
                raise ValueError(f"Restored tree did not match its verified inventory: {relative}")
            del records[relative]
            self._write_records(records)
            (self.root / record["archive"]).unlink()
        finally:
            shutil.rmtree(staging)

    def restore_all(self) -> list[str]:
        """Restore all archives; refuse modified existing paths or corruption."""
        records = self.records()
        restored = []
        for relative in list(records):
            self._restore_one(relative, records)
            restored.append(relative)
        return restored


@contextmanager
def output_lock(root: Path):
    root = root.resolve()
    directory = root / ".storage"
    if directory.is_symlink() or (directory / "run.lock").is_symlink():
        raise ValueError("Refusing a symlinked storage lock path")
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "run.lock").open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"Another pipeline/storage process holds the output lock: {root}") from error
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def space_usage(root):
    blocks = files = directories = 0
    seen = set()
    for directory, subdirs, names in os.walk(root, followlinks=False):
        directories += 1
        for name in subdirs + names:
            path = Path(directory) / name
            info = path.lstat()
            inode = (info.st_dev, info.st_ino)
            if inode not in seen:
                seen.add(inode)
                blocks += getattr(info, "st_blocks", (info.st_size + 511) // 512) * 512
            if name in names or path.is_symlink():
                files += 1
    return {"allocated_bytes": blocks, "files_and_symlinks": files, "directories": directories}

FINAL_MESSAGE = "TE library, Summary Figures, and TE Quantifications in Standard Formats Can Be Found in"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def layout(root, species=None):
    """Accept either an outer output directory or its species_EarlGrey child."""
    supplied = Path(root).expanduser().absolute()
    if supplied.is_symlink():
        raise ValueError(f"EarlGrey root must not use a symlinked directory: {supplied}")
    supplied = supplied.resolve()
    if not supplied.is_dir():
        raise ValueError(f"EarlGrey output directory does not exist: {supplied}")
    inferred = re.sub(r"_earlgrey$", "", supplied.name, flags=re.IGNORECASE)
    if species is None:
        if inferred == supplied.name:
            raise ValueError("Use --species when the root name does not end in _Earlgrey")
        species = inferred
    if not species or species in {".", ".."} or any(char in species for char in "/\\\n\r"):
        raise ValueError("Invalid EarlGrey species identifier")
    child = supplied / f"{species}_EarlGrey"
    is_inner = any((supplied / f"{species}_{suffix}").exists()
                   for suffix in ("Database", "summaryFiles", "RepeatModeler", "strainer"))
    if is_inner and child.exists():
        raise ValueError("Ambiguous nested EarlGrey output directories")
    if is_inner:
        inner, outer = supplied, supplied.parent
    else:
        inner, outer = child, supplied
    if inner.is_symlink():
        raise ValueError(f"EarlGrey inner directory must not be a symlink: {inner}")
    return outer, inner, species


def _required(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty completed EarlGrey output: {path}")


def _recorded_digest(store, relative):
    for record in store.records().values():
        for member in record["members"]:
            if member["name"] == relative and member["type"] == "file":
                return member["sha256"]
    return None


def _same_export(store, original, exported):
    _required(exported)
    actual = digest(original) if original.is_file() else _recorded_digest(
        store, original.relative_to(store.root).as_posix())
    if actual is None or actual != digest(exported):
        raise ValueError(f"EarlGrey summary export does not match its stage source: {exported}")


def validate_completed(inner: Path, species: str):
    """Require native completion and the exported products, retaining benchmarks."""
    if not inner.is_dir():
        raise ValueError(f"EarlGrey run directory does not exist: {inner}")
    store = ArchiveStore(inner)
    log = inner / f"{species}_EarlGrey.log"
    _required(log)
    text = log.read_text(errors="replace")
    completed_at = text.rfind(FINAL_MESSAGE)
    if completed_at < 0 or not re.search(r"<<<\s*Done in\s+", text[:completed_at]):
        raise ValueError("EarlGrey log has no complete final-stage completion receipt")
    trailing_start = text.find("\n", completed_at)
    trailing = text[trailing_start + 1:] if trailing_start >= 0 else ""
    if re.search(r"(?m)^\s*<<<\s*\S", trailing):
        raise ValueError("EarlGrey log contains a newer stage start after its completion receipt")
    if re.search(r"(?im)^\s*(?:ERROR\b|FATAL\b|Traceback\b|.*Error:)", trailing):
        raise ValueError("EarlGrey log reports a failure after its completion receipt")
    summary = inner / f"{species}_summaryFiles"
    for suffix in ("highLevelCount.txt", "highLevelCount.kable", "summaryPie.pdf",
                   "filteredRepeats.bed", "filteredRepeats.gff"):
        _required(summary / f"{species}.{suffix}")
    _required(summary / f"{species}_divergence_summary_table.tsv")
    baseline = inner / f"{species}_Database" / f"{species}-families.fa"
    _required(baseline)
    libraries = [inner / f"{species}_strainer" / f"{species}-families.fa{suffix}"
                 for suffix in (".strained", ".strained.clstrd.fa")]
    available = [path for path in libraries if path.is_file() and path.stat().st_size]
    if not available:
        raise ValueError("Missing completed top-level EarlGrey strained library")
    library = available[-1]
    _same_export(store, library, summary / library.name)
    masker = inner / f"{species}_RepeatMasker_Against_Custom_Library"
    pairs = [(path, path.with_suffix(".tbl")) for path in sorted(masker.glob("*.out"))
             if path.is_file() and path.with_suffix(".tbl").is_file()]
    if len(pairs) != 1:
        raise ValueError("Expected exactly one matching raw RepeatMasker .out/.tbl pair")
    for path in pairs[0]:
        _required(path)
    merged = inner / f"{species}_mergedRepeats"
    candidates = [merged / "looseMerge", merged]
    chosen = []
    for directory in candidates:
        try:
            for suffix in ("bed", "gff"):
                _same_export(store, directory / f"{species}.filteredRepeats.{suffix}",
                             summary / f"{species}.filteredRepeats.{suffix}")
        except ValueError:
            continue
        chosen.append(directory)
    if not chosen:
        raise ValueError("Final BED/GFF exports cannot be matched to loose or strict merge products")
    return {"log": str(log), "raw_repeatmasker_out": str(pairs[0][0]),
            "raw_repeatmasker_table": str(pairs[0][1]), "library": str(library),
            "summary_directory": str(summary), "merge_source": str(chosen[0])}


def plan_compaction(inner: Path, species: str):
    store = ArchiveStore(inner)
    records = store.records()
    candidates = []
    modeler = inner / f"{species}_RepeatModeler"
    candidates.extend(path for path in sorted(modeler.glob("RM_*")) if path.is_dir())
    strainer = inner / f"{species}_strainer"
    candidates.extend(path for path in sorted(strainer.glob("TS_*")) if path.is_dir())
    candidates.extend(inner / f"{species}_{name}" for name in ("mergedRepeats", "RepeatLandscape"))
    masker = inner / f"{species}_RepeatMasker_Against_Custom_Library"
    for pattern in ("*.align", "*.masked"):
        candidates.extend(sorted(masker.glob(pattern)))
    heliano = inner / f"{species}_heliano"
    for run in sorted(heliano.glob("HEL_*")):
        if not (run / "RC.representative.gff").is_file() or not (run / "RC.representative.gff").stat().st_size:
            continue
        candidates.extend(run / name for name in ("HLE1", "HLE2"))
        candidates.extend(path for path in sorted(run.glob("mmseqs*")) if path.is_dir())
    return [{"operation": "archive", "path": path.relative_to(inner).as_posix()}
            for path in candidates if (path.exists() or path.is_symlink())
            and path.relative_to(inner).as_posix() not in records]


def _preflight(inner, plan):
    """Check every candidate and inbound link before removing any loose files."""
    store = ArchiveStore(inner)
    targets = [inner / item["path"] for item in plan]
    for item in plan:
        store._inventory(item["path"])
    for directory, directories, files in os.walk(inner, followlinks=False):
        parent = Path(directory)
        directories[:] = [name for name in directories
                          if parent / name not in targets and name != ".storage"]
        for name in directories + files:
            path = parent / name
            if path.is_symlink():
                destination = path.resolve()
                if any(destination == target or target in destination.parents for target in targets):
                    raise ValueError(f"A retained symlink depends on an archive candidate: {path}")


def preview(inner, species):
    if not inner.is_dir():
        return {"read_only": True, "eligible": False, "errors": [f"EarlGrey run directory does not exist: {inner}"],
                "species": species, "outdir": str(inner), "plan": [], "existing_archives": [],
                "completion": None}
    store = ArchiveStore(inner)
    plan = plan_compaction(inner, species)
    errors = []
    completion = None
    try:
        completion = validate_completed(inner, species)
        for record in store.records().values():
            store._verify_archive(record)
        _preflight(inner, plan)
    except (OSError, ValueError, RuntimeError) as error:
        errors.append(str(error))
    return {"read_only": True, "eligible": not errors, "errors": errors,
            "species": species, "outdir": str(inner), "plan": plan,
            "existing_archives": list(store.records()), "completion": completion,
            "space": space_usage(inner)}


def compact(inner, species):
    """Caller owns the output lock and has established native-job inactivity."""
    completion = validate_completed(inner, species)
    destination = inner / ".storage/earlgrey_report.json"
    if destination.is_symlink():
        raise ValueError("Refusing symlinked EarlGrey storage report")
    plan = plan_compaction(inner, species)
    _preflight(inner, plan)
    store = ArchiveStore(inner)
    for record in store.records().values():
        store._verify_archive(record)
    before = space_usage(inner)
    actions = []
    for item in plan:
        record = store.archive(item["path"])
        actions.append({**item, "archive": record["archive"],
                        "members": len(record["members"])})
    report = {"read_only": False, "species": species, "outdir": str(inner),
              "completed_at_utc": datetime.now(timezone.utc).isoformat(),
              "completion": completion, "before": before, "after": space_usage(inner),
              "actions": actions, "archives": list(store.records())}
    destination = inner / ".storage/earlgrey_report.json"
    destination.parent.mkdir(exist_ok=True)
    if destination.is_symlink():
        raise ValueError("Refusing symlinked EarlGrey storage report")
    temporary = destination.with_suffix(f".{uuid.uuid4().hex}.partial")
    try:
        with temporary.open("x") as handle:
            json.dump(report, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return report


def _command_option(command, flag):
    found = [index for index, value in enumerate(command) if value == flag]
    if len(found) != 1 or found[0] + 1 >= len(command):
        raise ValueError(f"Wrapped EarlGrey command must specify {flag} exactly once")
    return command[found[0] + 1]


def _group_has_live_members(group):
    """Distinguish an unkillable zombie-only group on macOS from live work."""
    try:
        result = subprocess.run(["ps", "-axo", "pgid=,stat="], check=True,
                                capture_output=True, text=True)
    except (OSError, subprocess.SubprocessError):
        return None
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[0] == str(group) and not fields[1].startswith("Z"):
            return True
    return False


def _stop_native_group(process):
    # start_new_session owns a new group with the child's PID. Never signal
    # a caller's process group, including when using a PID-less test double.
    pid = getattr(process, "pid", None)
    if type(pid) is not int or pid <= 1 or pid == os.getpgrp():
        return

    def send(sig):
        try:
            os.killpg(pid, sig)
            return True
        except ProcessLookupError:
            return False
        except PermissionError:
            # macOS can return EPERM when every remaining member is a zombie.
            # Never dismiss a permission failure while live workers may exist.
            if process.poll() is not None and _group_has_live_members(pid) is False:
                return False
            raise

    send(signal.SIGTERM)
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        send(signal.SIGKILL)
        process.wait(timeout=10)
    else:
        # A leader can exit before its children. Signal our known owned group
        # directly; signal-0 probes can be denied for reaped/zombie groups.
        send(signal.SIGKILL)


def run_native(command):
    """Keep native descendants inside our lock through controlled shutdown."""
    process = None
    install_handler = threading.current_thread() is threading.main_thread()
    previous_handler = signal.getsignal(signal.SIGTERM) if install_handler else None

    def terminated(signum, frame):
        raise SystemExit(128 + signum)

    if install_handler:
        signal.signal(signal.SIGTERM, terminated)
    try:
        process = subprocess.Popen(command, start_new_session=True)
        returncode = process.wait()
        if returncode:
            _stop_native_group(process)
        return returncode
    except BaseException:
        if install_handler:
            signal.signal(signal.SIGTERM, signal.SIG_IGN)
        if process is not None:
            try:
                _stop_native_group(process)
            except BaseException as error:
                logging.getLogger(__name__).warning("Could not fully stop EarlGrey workers: %s", error)
                pid = getattr(process, "pid", None)
                if type(pid) is int and pid > 1 and pid != os.getpgrp():
                    try:
                        os.killpg(pid, signal.SIGKILL)
                        process.wait(timeout=10)
                    except BaseException as cleanup_error:
                        logging.getLogger(__name__).warning("Forced EarlGrey cleanup failed: %s", cleanup_error)
        raise
    finally:
        if install_handler:
            signal.signal(signal.SIGTERM, previous_handler)


def run_wrapped(outer, inner, species, command):
    if command and command[0] == "--":
        command = command[1:]
    if not command or Path(command[0]).name.lower() != "earlgrey":
        raise ValueError("--run must be followed by an explicit earlGrey command")
    if any(value.startswith(("-s", "-o", "-g")) and value not in {"-s", "-o", "-g"}
           for value in command[1:]):
        raise ValueError("Use separate -g, -s and -o values; combined option forms cannot be validated")
    if _command_option(command, "-s") != species:
        raise ValueError("Wrapped EarlGrey -s differs from the selected species")
    output = Path(_command_option(command, "-o")).expanduser().resolve()
    if output != outer:
        raise ValueError(f"Wrapped EarlGrey -o must be the outer directory: {outer}")
    _command_option(command, "-g")
    with output_lock(inner):
        ArchiveStore(inner).restore_all()
        returncode = run_native(command)
        if returncode:
            return {"returncode": returncode, "compacted": False,
                    "reason": "EarlGrey did not exit successfully; working files retained"}
        return {"returncode": 0, "compacted": True, "storage": compact(inner, species)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", required=True, type=Path, help="Outer species_Earlgrey or inner species_EarlGrey directory")
    parser.add_argument("--species", help="Infer from root name when omitted")
    parser.add_argument("--apply", action="store_true", help="Apply compaction/restore to an idle run")
    parser.add_argument("--idle", action="store_true", help="Confirm no native EarlGrey jobs/workers use this output")
    parser.add_argument("--restore", action="store_true")
    parser.add_argument("--run", nargs=argparse.REMAINDER, help="Restore, run the following earlGrey command, then compact successful results")
    args = parser.parse_args(argv)
    try:
        if args.run is not None and (args.apply or args.restore):
            raise ValueError("Use --run separately from --apply/--restore")
        if args.apply and not args.idle:
            raise ValueError("--apply requires --idle: native EarlGrey jobs do not hold this storage lock")
        # New wrapped runs may start with an empty, explicit outer output root.
        if args.run is not None and not args.root.exists():
            if args.root.is_symlink():
                raise ValueError("Refusing symlinked output root")
            args.root.mkdir(parents=True)
        outer, inner, species = layout(args.root, args.species)
        if args.run is not None:
            if inner.is_dir() and any(inner.iterdir()) and not args.idle:
                raise ValueError("--run on an existing run requires --idle: native EarlGrey jobs do not hold this storage lock")
            result = run_wrapped(outer, inner, species, args.run)
        elif args.restore:
            store = ArchiveStore(inner)
            if args.apply:
                with output_lock(inner):
                    result = {"restored": store.restore_all(), "read_only": False}
            else:
                result = {"restore": list(store.records()), "read_only": True, "outdir": str(inner)}
        elif args.apply:
            with output_lock(inner):
                result = compact(inner, species)
        else:
            result = preview(inner, species)
        print(json.dumps(result, indent=2))
        return int(result.get("returncode", 0))
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(2, f"{parser.prog}: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
