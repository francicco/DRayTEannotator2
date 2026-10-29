import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence, Optional


def run_command_to_file(
    cmd: Sequence[str],
    log_file: str | Path,
    logger,
    cwd: str | Path | None = None,
    env: Optional[dict] = None,
    label: str = "Command",
) -> None:
    """Write both output streams directly to disk, with brief status messages."""
    log_file = Path(log_file).resolve()
    log_file.parent.mkdir(parents=True, exist_ok=True)
    command = [str(value) for value in cmd]
    cmd_str = shlex.join(command)
    logger.info("Running %s; full output: %s", label, log_file)

    # Keep earlier attempts available when a failed or outdated stage is rerun.
    with log_file.open("a", encoding="utf-8") as logh:
        timestamp = datetime.now(timezone.utc).isoformat()
        logh.write(f"\n=== {label} started {timestamp} ===\nCommand: {cmd_str}\n")
        logh.flush()
        try:
            subprocess.run(
                command,
                check=True,
                cwd=str(cwd) if cwd is not None else None,
                env=env,
                stdout=logh,
                stderr=subprocess.STDOUT,
            )
        except (subprocess.CalledProcessError, OSError) as exc:
            if isinstance(exc, subprocess.CalledProcessError):
                reason = f"exit code {exc.returncode}"
            else:
                reason = str(exc)
            message = f"{label} failed ({reason}); see full log: {log_file}"
            logh.write(f"\n=== {message} ===\n")
            logger.error("%s", message)
            raise RuntimeError(message) from exc
        logh.write(f"\n=== {label} completed successfully ===\n")

    logger.info("%s completed successfully; full output: %s", label, log_file)


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
