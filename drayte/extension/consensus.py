from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


LOGGER = logging.getLogger(__name__)


@dataclass
class ExtensionResult:
    te_id: str
    rep_fa: Optional[Path]
    msa_fa: Optional[Path]
    png_file: Optional[Path]
    hit_count: int
    consensus_length: int
    category: str


def run_extend_consensus(
    extend_script: str | Path,
    genome_2bit: str | Path,
    family_fasta: str | Path,
    workdir: str | Path,
    logfile: str | Path,
) -> None:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "perl",
        str(extend_script),
        "-genome",
        str(genome_2bit),
        "-family",
        str(family_fasta),
        "-outdir",
        ".",
    ]

    LOGGER.info("Running extension: %s", " ".join(cmd))
    with open(logfile, "w") as logh:
        process = subprocess.Popen(
            cmd,
            cwd=str(workdir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            line = line.rstrip()
            if line:
                LOGGER.info("[extend] %s", line)
                logh.write(line + "\n")
        rc = process.wait()
        if rc != 0:
            raise RuntimeError(f"Extension tool failed for {family_fasta} with exit code {rc}")


def postprocess_extension_outputs(
    te_id: str,
    te_workdir: str | Path,
    ) -> tuple[Optional[Path], Optional[Path], Optional[Path], int, int]:
    
    te_workdir = Path(te_workdir)

    rep_src = te_workdir / "rep"
    msa_src = te_workdir / "MSA-extended_with_rmod_cons.fa"
    png_src = te_workdir / "img.png"
    unextended = te_workdir / "repseq.unextended"

    rep_out = te_workdir / f"{te_id}_rep.fa"
    msa_out = te_workdir / f"{te_id}_MSA_extended.fa"
    png_out = te_workdir / f"{te_id}.png"

    hit_count = 0
    if unextended.exists():
        with open(unextended) as handle:
            hit_count = sum(
                1 for line in handle
                if line.startswith(">")
                and not line.startswith(">CONSENSUS-")
                and not line.startswith(">CORECONS")
                and not line.startswith(">repam-newrep")
            )

    consensus_length = 0

    if rep_src.exists():
        text = rep_src.read_text()
        text = text.replace("repam-newrep", te_id)
        rep_out.write_text(text)

        with open(rep_out) as handle:
            consensus_length = sum(len(line.strip()) for line in handle if not line.startswith(">"))

    if msa_src.exists():
        text = msa_src.read_text()
        text = text.replace("repam-newrep", te_id)
        text = text.replace("CORECONS", f"CONSENSUS-{te_id}")
        msa_out.write_text(text)

    if png_src.exists():
        shutil.copy2(png_src, png_out)

    return (
        rep_out if rep_out.exists() else None,
        msa_out if msa_out.exists() else None,
        png_out if png_out.exists() else None,
        hit_count,
        consensus_length,
    )


def categorize_extension(hit_count: int, consensus_length: int) -> str:
    if hit_count < 10:
        return "reject"
    if consensus_length > 15000:
        return "possible_SD"
    return "likely_TE"
