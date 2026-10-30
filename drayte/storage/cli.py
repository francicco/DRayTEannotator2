"""Preview, compact, or restore one idle DrayTE output tree."""
import argparse
import json
import logging

from drayte.pipeline.config import load_config
from .archive import ArchiveStore
from .policy import compact_run, output_lock, plan_compaction, validate_completed_run


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Original pipeline YAML config")
    parser.add_argument("--restore", action="store_true", help="Restore managed archives instead of compacting")
    parser.add_argument("--apply", action="store_true", help="Apply the plan; only use on an idle output directory")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    root = config.outdir_path.resolve()
    if not root.is_dir():
        parser.error(f"Output directory does not exist: {root}")
    store = ArchiveStore(root)
    if not args.apply:
        plan = list(store.records()) if args.restore else plan_compaction(config)
        validation = "restore preview" if args.restore else "already compacted" if store.records() else "not checked"
        if not args.restore and not store.records():
            try:
                validate_completed_run(config)
                validation = "completed outputs validated"
            except (OSError, ValueError, KeyError) as error:
                validation = f"blocked: {error}"
        print(json.dumps({"mode": "restore" if args.restore else "compact", "apply": False,
                          "outdir": str(root), "plan": plan,
                          "validation": validation, "existing_archives": list(store.records())}, indent=2))
        print("Preview only. No output files changed. Validate job inactivity before --apply.")
        return 0
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    logger = logging.getLogger("drayte.storage")
    with output_lock(root):
        if args.restore:
            restored = store.restore_all()
            logger.info("Restored %s managed paths", len(restored))
        else:
            compact_run(config, logger)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
