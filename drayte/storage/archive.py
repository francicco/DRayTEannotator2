"""Verified, recoverable archives for pipeline intermediates.

The caller must hold the run's storage lock. Archive records are committed before
source removal. Interrupted removal is recoverable by merging missing files from
the archive; changed existing files are never overwritten.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import stat
import tarfile
import uuid
from pathlib import Path, PurePosixPath


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

    def archive(self, relative: str | Path) -> dict:
        """Verify a compressed archive, commit its record, then remove source."""
        relative = self._relative(relative)
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
        if source.is_dir():
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
