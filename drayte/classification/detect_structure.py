from __future__ import annotations

import argparse

from .structure_detect import (
    detect_tirs_from_fasta,
    write_tir_structure_tsv,
)


def main():

    parser = argparse.ArgumentParser(
        description="Detect TE structural signals"
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Consensus FASTA"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Structure evidence TSV"
    )

    parser.add_argument(
        "--window",
        type=int,
        default=250,
    )

    parser.add_argument(
        "--min-len",
        type=int,
        default=15,
    )

    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.80,
    )

    args = parser.parse_args()

    detections = detect_tirs_from_fasta(
        fasta=args.fasta,
        window=args.window,
        min_len=args.min_len,
        min_identity=args.min_identity,
    )

    write_tir_structure_tsv(
        detections,
        args.output,
    )

    n = sum(d.tir_present for d in detections)

    print(f"Detected TIRs in {n} families")


if __name__ == "__main__":
    main()
