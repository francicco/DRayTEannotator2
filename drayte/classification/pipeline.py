from __future__ import annotations

from pathlib import Path

from .hmmer import (
    parse_domtblout,
    run_hmmscan,
)
from .orfs import write_translated_orfs


def run_domain_annotation(
    consensus_fasta: str | Path,
    hmm_db: str | Path,
    outdir: str | Path,
    hmmscan_bin: str = "hmmscan",
    min_orf_nt: int = 500,
    include_reverse_orfs: bool = True,
    cpu: int = 1,
):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    orf_fasta = outdir / "translated_orfs.fa"
    domtblout = outdir / "domains.domtblout"

    #
    # ORFs
    #

    write_translated_orfs(
        input_fasta=consensus_fasta,
        output_fasta=orf_fasta,
        minsize_nt=min_orf_nt,
        include_reverse=include_reverse_orfs,
    )

    #
    # hmmscan
    #

    run_hmmscan(
        hmm_db=hmm_db,
        proteins_fasta=orf_fasta,
        domtblout=domtblout,
        hmmscan_bin=hmmscan_bin,
        cpu=cpu,
    )

    #
    # parse
    #

    return parse_domtblout(domtblout)
