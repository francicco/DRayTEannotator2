"""Storage policy for completed pipeline runs; never prune unfinished evidence."""
from __future__ import annotations

import fcntl
import json
import os
import shutil
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from .archive import ArchiveStore


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


def storage_mode(config):
    settings = config.extra.get("storage")
    settings = {} if settings is None else settings
    if not isinstance(settings, dict):
        raise ValueError("storage must be a YAML mapping")
    mode = settings.get("mode", "keep")
    if mode not in {"keep", "compact"}:
        raise ValueError("storage.mode must be keep or compact")
    return mode


def read_manifest(path):
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Invalid manifest: {path}")
    return value


def require_files(paths):
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Storage compaction requires a complete run; missing/empty: {path}")


def repeatmasker_contexts(config):
    root = config.outdir_path.resolve()
    species = config.species
    executable = str(config.extra.get("repeatmasker_bin", "RepeatMasker"))
    final = read_manifest(root / "classification/classification.manifest.json")
    for stage, relative, library, flags in (
        ("initial", "discovery/rmasker_dir", root / "discovery/rmodeler_dir/consensi.fa.classified",
         ["-xsmall", "-gc", "-gff", "-norna", "-s"]),
        ("final", "final_annotation", Path(final["classified_library"]),
         ["-s", "-xsmall", "-gc", "-gff", "-norna", "-a"]),
    ):
        directory = root / relative
        suffixes = (".out", ".tbl") if stage == "initial" else (".out", ".tbl", ".masked", ".cat.gz", ".align")
        yield stage, directory, library, [directory / f"{species}.fa{suffix}" for suffix in suffixes], [executable, *flags]


def adopt_legacy_repeatmasker_records(config):
    """Persist explicitly labelled legacy adoption before archiving outputs."""
    from drayte.utils.repeatmasker_checkpoint import input_signature, write_checkpoint
    root = config.outdir_path.resolve()
    for _, directory, library, outputs, options in repeatmasker_contexts(config):
        marker = directory / "repeatmasker.checkpoint.json"
        if not marker.exists():
            write_checkpoint(marker, input_signature(root / f"discovery/assemblies_dir/{config.species}.fa",
                                                    library, options), outputs, legacy=True)


def validate_completed_run(config):
    """Validate persisted stage outputs and classification biological digests.

    This is not a substitute for the pipeline's config-specific resume checks.
    It proves the products to be archived still match the completed run.
    """
    from drayte.pipeline import classification, extension, reclassify

    root = config.outdir_path.resolve()
    species = config.species
    require_files([
        root / "discovery/rmodeler_dir/consensi.fa.classified",
        root / f"discovery/rmasker_dir/{species}.fa.out",
        root / f"discovery/rmasker_dir/{species}.fa.tbl",
        root / f"extension/blastfiles/{species}-families.mod.fa",
        root / f"reclassify/{species}_extended_rep.fa",
        root / f"reclassify/{species}_extended_rep.fa.classified",
        root / "curation/Final.RepeatModeler.Lib.fa",
        *[root / f"final_annotation/{species}.fa{suffix}" for suffix in
          (".out", ".tbl", ".masked", ".cat.gz", ".align")],
        root / f"Defragmenting_RModAnnotation/{species}.annotation_refinement.tsv",
        root / f"Defragmenting_RModAnnotation/{species}.filteredRepeats.tsv",
        root / f"Defragmenting_RModAnnotation/{species}.annotation_refinement.stats.tsv",
        root / f"{species}.Defragmenting_RModAnnotation_noFiltered/{species}.highLevelCount.txt",
        root / f"{species}.Defragmenting_RModAnnotation_Filtered/{species}.highLevelCount.txt",
    ])
    result = read_manifest(root / "extension/extension.manifest.json")
    valid, reason = extension.validate_extension_outputs(
        result, root / "extension/extensionwork", root / "extension/final_consensuses"
    )
    if not valid:
        raise ValueError(f"Cannot compact an incomplete extension: {reason}")
    checkpoint = root / "reclassify/repeatclassifier.checkpoint.json"
    if not reclassify.checkpoint_is_valid(
        checkpoint, root / f"reclassify/{species}_extended_rep.fa",
        root / f"reclassify/{species}_extended_rep.fa.classified",
        config.extra.get("repeatclassifier_bin", "RepeatClassifier"),
    ):
        raise ValueError("Cannot compact: RepeatClassifier completion checkpoint is missing or invalid")
    result = read_manifest(root / "classification/classification.manifest.json")
    if not classification.manifest_outputs_are_valid(result):
        raise ValueError("Cannot compact: classification outputs are incomplete")
    signature = result.get("input_signature") or {}
    input_library = Path(result.get("input_library", root / "curation/Final.RepeatModeler.Lib.fa"))
    structural = root / f"classification/{species}.structure_input.fa"
    expected = {
        "library_sha256": classification.sha256_file(input_library),
        "structure_library_sha256": classification.sha256_file(structural),
        "extensionwork_sha256": classification.selected_extensionwork_sha256(root / "extension/extensionwork"),
    }
    if any(signature.get(key) != value for key, value in expected.items()):
        raise ValueError("Cannot compact: classification biological fingerprints do not match")
    if not classification.classified_library_matches_input(structural, Path(result["classified_library"])):
        raise ValueError("Cannot compact: classified family sequences do not match their input")
    required_manifest_fields = {
        "discovery": "raw_library", "curation": "final_library",
        "final_annotation": "repeatmasker_out",
    }
    for stage, key in required_manifest_fields.items():
        data = read_manifest(root / stage / f"{stage}.manifest.json")
        if not isinstance(data.get(key), str):
            raise ValueError(f"Incomplete {stage} manifest: missing {key}")
        require_files([Path(data[key])])
    refinement = read_manifest(root / f"Defragmenting_RModAnnotation/{species}.annotation_refinement.manifest.json")
    if "refinement_algorithm_version" not in refinement:
        raise ValueError("Cannot compact: refinement completion manifest is incomplete")
    inspection = config.extra.get("family_inspection") or {}
    if inspection.get("enabled", True):
        manifest = read_manifest(root / "family_inspection/family_inspection.manifest.json")
        require_files([Path(manifest["manifest"])])
    # Check modern completion records and their inputs/outputs; a state label
    # alone is not proof of completion. Legacy adoption remains explicit.
    from drayte.utils.repeatmasker_checkpoint import checkpoint_is_valid, input_signature, legacy_outputs_usable
    for stage, directory, library, outputs, options in repeatmasker_contexts(config):
        marker = directory / "repeatmasker.checkpoint.json"
        if marker.exists():
            signature = input_signature(root / f"discovery/assemblies_dir/{species}.fa",
                                        library, options)
            if not checkpoint_is_valid(marker, signature, outputs):
                raise ValueError(f"Cannot compact an invalid RepeatMasker completion record: {directory}")
        elif not legacy_outputs_usable(outputs, directory / f"{species}.RepeatMasker.log"):
            raise ValueError(f"Cannot compact incomplete legacy RepeatMasker outputs: {directory}")
        elif stage == "final" and outputs[0].stat().st_mtime_ns < library.stat().st_mtime_ns:
            raise ValueError("Cannot compact legacy final annotation: classified library is newer than its hits")


def plan_compaction(config, *, runtime_inputs_known=False):
    root = config.outdir_path.resolve()
    species = config.species
    records = ArchiveStore(root).records()
    archive = [
        "extension/extensionwork", "extension/extract_align", "extension/extendlogs",
        "extension/images_and_alignments", "curation/te-aid",
        "classification/pfam/chunks", "classification/pfam/hmmscan_chunks",
        "classification/dfam/chunks", "classification/dfam/tblouts",
        f"discovery/rmasker_dir/{species}.fa.masked",
        f"final_annotation/{species}.fa.masked", f"final_annotation/{species}.fa.align",
    ]
    for category in ("DNA", "LINE", "SINE", "LTR", "RC", "Penelope", "Unclassified", "Unknown", "Other"):
        archive.append(f"family_inspection/{category}")
    plan = [{"operation": "archive", "path": rel}
            for rel in archive if rel not in records and (root / rel).exists()]
    modeler = root / "discovery/rmodeler_dir"
    if modeler.is_dir():
        plan.extend({"operation": "delete_scratch", "path": str(path.relative_to(root))}
                    for path in sorted(modeler.glob("RM_*")) if path.is_dir())
    plan.extend({"operation": "delete_scratch", "path": rel}
                for rel in ("classification/mmseqs/tmp", "heliano/mmseqs_tmp")
                if (root / rel).is_dir())
    if not runtime_inputs_known:
        # An old YAML need not contain the CLI-only library inputs originally
        # used. Preserve scratch reversibly rather than deleting a possible
        # source required by the original resume command.
        for item in plan:
            if item["operation"] == "delete_scratch":
                item["operation"] = "archive"
                item["reason"] = "Original runtime-only inputs are unknown; preserve reversibly"
    return plan


def cleanup_stage_scratch(config, stage, logger):
    """Reduce peak storage only after the producing stage has succeeded.

    Archive operations wait until full completion, since downstream stages still
    need loose evidence. This function never removes evidence or checkpoints.
    """
    if storage_mode(config) != "compact":
        return []
    root = config.outdir_path.resolve()
    species = config.species
    if stage == "discovery":
        require_files([root / "discovery/rmodeler_dir/consensi.fa.classified",
                       root / f"discovery/rmasker_dir/{species}.fa.out",
                       root / f"discovery/rmasker_dir/{species}.fa.tbl"])
        candidates = sorted((root / "discovery/rmodeler_dir").glob("RM_*"))
    elif stage == "heliano":
        if not (root / f"heliano/{species}.heliano.unique_vs_curated.fa").is_file():
            return []
        require_files([root / f"heliano/{species}.heliano.fa"])
        candidates = [root / "heliano/mmseqs_tmp"]
    elif stage == "classification":
        from drayte.pipeline.classification import manifest_outputs_are_valid
        manifest = read_manifest(root / "classification/classification.manifest.json")
        if not manifest_outputs_are_valid(manifest):
            raise ValueError("Cannot prune scratch for incomplete classification")
        candidates = [root / "classification/mmseqs/tmp"]
    else:
        raise ValueError(f"Unsupported scratch-cleanup stage: {stage}")
    targets = [path for path in candidates if path.is_dir() or path.is_symlink()]
    assert_safe_deletions(config, targets)
    if targets and not shutil.rmtree.avoids_symlink_attacks:
        raise ValueError("This Python does not support symlink-resistant scratch deletion")
    removed = []
    journal = root / ".storage/scratch_actions.jsonl"
    journal.parent.mkdir(parents=True, exist_ok=True)
    if journal.is_symlink():
        raise ValueError("Refusing a symlinked scratch journal")
    for path in targets:
        action = {"stage": stage, "path": str(path.relative_to(root)),
                  "at_utc": datetime.now(timezone.utc).isoformat()}
        with journal.open("a") as handle:
            handle.write(json.dumps({**action, "state": "planned"}) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        shutil.rmtree(path)
        with journal.open("a") as handle:
            handle.write(json.dumps({**action, "state": "removed"}) + "\n")
        logger.info("Storage: removed completed %s scratch %s", stage, path)
        removed.append(str(path.relative_to(root)))
    return removed


def assert_safe_deletions(config, targets):
    """Reject links or configured external inputs that depend on scratch."""
    root = config.outdir_path.resolve()
    for target in targets:
        for path in (target, *target.parents):
            if path == root:
                break
            if path.is_symlink():
                raise ValueError(f"Refusing symlinked scratch path: {path}")
    def protected(path):
        resolved = path.resolve()
        return any(resolved == target or target in resolved.parents for target in targets)
    def visit_settings(value):
        if isinstance(value, dict):
            for child in value.values():
                visit_settings(child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                visit_settings(child)
        elif isinstance(value, str) and value and protected(Path(value).expanduser()):
            raise ValueError(f"A configured input references scratch selected for deletion: {value}")
    visit_settings(config.extra)
    visit_settings(str(config.genome_path))
    visit_settings(str(root / "discovery/rmodeler_dir/consensi.fa.classified"))
    def walk_error(error):
        raise error
    for directory, directories, files in os.walk(root, followlinks=False, onerror=walk_error):
        parent = Path(directory)
        directories[:] = [name for name in directories
                          if parent / name not in targets and name != ".storage"]
        for name in directories + files:
            path = parent / name
            if path.is_symlink() and protected(path):
                raise ValueError(f"A retained symlink depends on scratch: {path}")


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


def compact_run(config, logger, *, runtime_inputs_known=False):
    """Caller holds output_lock. Default CLI use is a read-only plan."""
    root = config.outdir_path.resolve()
    plan = plan_compaction(config, runtime_inputs_known=runtime_inputs_known)
    if not plan:
        logger.info("Storage: no new compaction candidates")
        return {"actions": [], "archives": ArchiveStore(root).records()}
    if ArchiveStore(root).records():
        raise ValueError("Restore existing managed archives before compacting a partially compacted run")
    validate_completed_run(config)
    targets = [root / item["path"] for item in plan if item["operation"] == "delete_scratch"]
    assert_safe_deletions(config, targets)
    if targets and not shutil.rmtree.avoids_symlink_attacks:
        raise ValueError("This Python does not support symlink-resistant scratch deletion")
    before = space_usage(root)
    adopt_legacy_repeatmasker_records(config)
    actions = []
    store = ArchiveStore(root)
    # Free proven-disposable scratch first, making headroom for archive creation.
    for target in targets:
        shutil.rmtree(target)
        actions.append({"operation": "delete_scratch", "path": str(target.relative_to(root))})
        logger.info("Storage: removed completed-stage scratch %s", target)
    for item in plan:
        if item["operation"] == "archive":
            record = store.archive(item["path"])
            actions.append({**item, "archive": record["archive"],
                            "archive_sha256": record["archive_sha256"], "members": len(record["members"])})
            logger.info("Storage: archived %s", item["path"])
    report = {
        "species": config.species,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "before": before, "after": space_usage(root),
        "actions": actions,
        "archives": {path: {"archive": record["archive"], "archive_sha256": record["archive_sha256"],
                            "members": len(record["members"])} for path, record in store.records().items()},
        "notes": ["Restore before opening loose per-family evidence paths.",
                  "Only completed-stage scratch was deleted; evidence is archived."]
    }
    report_path = root / ".storage/storage_report.json"
    temporary = report_path.with_suffix(f".{uuid.uuid4().hex}.partial")
    with temporary.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, report_path)
    logger.info("Storage: compacted files %s -> %s; allocated bytes %s -> %s",
                before["files_and_symlinks"], report["after"]["files_and_symlinks"],
                before["allocated_bytes"], report["after"]["allocated_bytes"])
    return report
