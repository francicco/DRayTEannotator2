import argparse

from .classify import classify_family
from .features import build_families_from_evidence
from .hmmer import parse_domtblout
from .io import load_families_tsv, write_classification_tsv


def classify_families(families):

    results = []

    for f in families:

        c = classify_family(f)

        results.append({
            "family_id": f.family_id,
            **c
        })

    return results


def main():

    parser = argparse.ArgumentParser(
        description="DRayTE evidence-based TE classifier"
    )

    input_group = parser.add_mutually_exclusive_group(required=True)

    input_group.add_argument(
        "--input",
        help="Input family TSV"
    )

    input_group.add_argument(
        "--fasta",
        help="Consensus FASTA"
    )

    parser.add_argument(
        "--domtblout",
        help="hmmscan domtblout file",
        default=None,
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output classification TSV"
    )

    parser.add_argument(
        "--min-orf-nt",
        type=int,
        default=500,
        help="Minimum ORF size in nucleotides"
    )

    parser.add_argument(
        "--forward-only",
        action="store_true",
        help="Do not scan reverse strand ORFs"
    )

    args = parser.parse_args()

    #
    # TSV mode
    #

    if args.input:

        families = load_families_tsv(args.input)

    #
    # FASTA mode
    #

    else:

        if args.domtblout:
            domain_hits = parse_domtblout(args.domtblout)
        else:
            domain_hits = []

        families = build_families_from_evidence(
            consensus_fasta=args.fasta,
            domain_hits=domain_hits,
            min_orf_nt=args.min_orf_nt,
            include_reverse_orfs=not args.forward_only,
        )

    results = classify_families(families)

    write_classification_tsv(
        results,
        args.output
    )

    print(f"Classified {len(results)} families")


if __name__ == "__main__":
    main()
