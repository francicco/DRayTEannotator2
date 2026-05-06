from __future__ import annotations

import argparse
import csv


def load_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def evaluate(expected_rows, predicted_rows):

    pred_by_id = {
        row["family_id"]: row
        for row in predicted_rows
    }

    total = 0

    class_ok = 0
    order_ok = 0
    superfamily_ok = 0

    missing = 0

    for exp in expected_rows:

        family_id = exp["family_id"]

        if family_id not in pred_by_id:
            missing += 1
            continue

        pred = pred_by_id[family_id]

        total += 1

        if pred["class"] == exp["expected_class"]:
            class_ok += 1

        if pred["order"] == exp["expected_order"]:
            order_ok += 1

        if pred["superfamily"] == exp["expected_superfamily"]:
            superfamily_ok += 1

    return {
        "total_compared": total,
        "missing_predictions": missing,
        "class_accuracy": class_ok / total if total else 0,
        "order_accuracy": order_ok / total if total else 0,
        "superfamily_accuracy": superfamily_ok / total if total else 0,
    }


def main():

    parser = argparse.ArgumentParser(
        description="Evaluate DRayTE TE classifications"
    )

    parser.add_argument(
        "--classifications",
        required=True,
    )

    parser.add_argument(
        "--expected",
        required=True,
    )

    args = parser.parse_args()

    expected = load_tsv(args.expected)
    predicted = load_tsv(args.classifications)

    metrics = evaluate(expected, predicted)

    print()
    print("=== DRayTE Classification Evaluation ===")
    print()

    for k, v in metrics.items():

        if isinstance(v, float):
            print(f"{k}: {v:.4f}")
        else:
            print(f"{k}: {v}")

    print()


if __name__ == "__main__":
    main()
