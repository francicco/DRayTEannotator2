from pathlib import Path
import argparse

from drayte.reporting.SummaryFilesGen import run_summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rmout", required=True, type=Path)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--outdir", required=True, type=Path)

    args = parser.parse_args()

    genome_size = 0
    with open(args.genome) as f:
        for line in f:
            if not line.startswith(">"):
                genome_size += len(line.strip())

    run_summary(
        rmout=args.rmout,
        genome_size=genome_size,
        species=args.species,
        outdir=args.outdir,
    )


if __name__ == "__main__":
    main()
