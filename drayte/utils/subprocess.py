import shlex
import subprocess
from pathlib import Path
from typing import Sequence


def run_command(
    cmd: Sequence[str],
    logger,
    cwd: str | Path | None = None,
    env: dict | None = None,
) -> None:
    cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)
    logger.info("Running command: %s", cmd_str)

    result = subprocess.run(
        [str(x) for x in cmd],
        cwd=str(cwd) if cwd else None,
        env=env,
        capture_output=True,
        text=True,
    )

    if result.stdout:
        logger.info(result.stdout.rstrip())

    if result.returncode != 0:
        if result.stderr:
            logger.error(result.stderr.rstrip())
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {cmd_str}")
