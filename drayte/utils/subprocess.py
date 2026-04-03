import shlex
import subprocess
from pathlib import Path
from typing import Sequence, Optional


def run_command(
    cmd: Sequence[str],
    logger,
    cwd: str | Path | None = None,
    env: Optional[dict] = None,
    prefix: Optional[str] = None,
) -> None:
    cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)
    logger.info("Running command: %s", cmd_str)

    process = subprocess.Popen(
        [str(x) for x in cmd],
        cwd=str(cwd) if cwd else None,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    assert process.stdout is not None

    for line in process.stdout:
        line = line.rstrip()
        if not line:
            continue
        if prefix:
            logger.info("[%s] %s", prefix, line)
        else:
            logger.info(line)

    returncode = process.wait()
    if returncode != 0:
        raise RuntimeError(f"Command failed with exit code {returncode}: {cmd_str}")
