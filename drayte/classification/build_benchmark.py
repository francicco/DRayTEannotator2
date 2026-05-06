from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO

from .parsers.repeatmasker import (
    parse_repeatmasker_header,
    classification_from_repeatmasker,
)


def build_expected_from_fasta(fasta: str | Path, output: str | Path) -> int:
    rows = []

    for rec in SeqIO.parse(str(fasta), "fasta"):
        parsed = parse_repeatmasker_header(rec.id)

        if parsed is None:
            continue

        rows.append(classification_from_repeatmasker(parsed))

    with open(output, "w") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=[
                "family_id",
                "expected_class",
                "expected_order",
                "expected_superfamily",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    return len(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Build DRayTE classification benchmark truth table from curated RepeatMasker-style FASTA"
    )
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    n = build_expected_from_fasta(args.fasta, args.output)
    print(f"Wrote {n} expected classifications")


if __name__ == "__main__":
    main()
