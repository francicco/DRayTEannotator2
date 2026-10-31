"""Preview and reversibly compact idle EarlGrey output, or wrap EarlGrey runs.

Native EarlGrey does not know about these archives. Restore before invoking it
again, or invoke it through ``--run`` so restoration happens automatically.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import signal
import subprocess
import threading
import uuid
from datetime import datetime, timezone
from pathlib import Path

from .archive import ArchiveStore
from .policy import output_lock, space_usage

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
    parser = argparse.ArgumentParser(description=__doc__)
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
        parser.exit(2, f"earlgrey-storage: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
