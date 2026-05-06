import argparse

from .io import load_families_tsv, write_classification_tsv
from .classify import classify_family


def main():

    parser = argparse.ArgumentParser(
        description="DRayTE evidence-based TE classifier"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input family TSV"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output classification TSV"
    )

    args = parser.parse_args()

    families = load_families_tsv(args.input)

    results = []

    for f in families:

        c = classify_family(f)

        results.append({
            "family_id": f.family_id,
            **c
        })

    write_classification_tsv(
        results,
        args.output
    )

    print(f"Classified {len(results)} families")


if __name__ == "__main__":
    main()
