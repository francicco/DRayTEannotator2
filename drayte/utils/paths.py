from pathlib import Path


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def stage_dir(outdir: str | Path, stage_name: str) -> Path:
    return ensure_dir(Path(outdir) / stage_name)
