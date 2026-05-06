from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict


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

    class_confusion = Counter()
    order_confusion = Counter()

    for exp in expected_rows:

        family_id = exp["family_id"]

        if family_id not in pred_by_id:
            missing += 1
            continue

        pred = pred_by_id[family_id]

        total += 1

        exp_class = exp["expected_class"]
        pred_class = pred["class"]

        exp_order = exp["expected_order"]
        pred_order = pred["order"]

        class_confusion[(exp_class, pred_class)] += 1
        order_confusion[(exp_order, pred_order)] += 1

        if pred_class == exp_class:
            class_ok += 1

        if pred_order == exp_order:
            order_ok += 1

        if pred["superfamily"] == exp["expected_superfamily"]:
            superfamily_ok += 1

    return {
        "metrics": {
            "total_compared": total,
            "missing_predictions": missing,
            "class_accuracy": class_ok / total if total else 0,
            "order_accuracy": order_ok / total if total else 0,
            "superfamily_accuracy": superfamily_ok / total if total else 0,
        },
        "class_confusion": class_confusion,
        "order_confusion": order_confusion,
    }


def print_confusion(title, confusion, top_n=25):

    print()
    print(f"=== {title} ===")
    print()

    for (expected, predicted), n in confusion.most_common(top_n):
        print(f"{n:5d}  expected={expected:<15} predicted={predicted}")


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

    results = evaluate(expected, predicted)

    metrics = results["metrics"]

    print()
    print("=== DRayTE Classification Evaluation ===")
    print()

    for k, v in metrics.items():

        if isinstance(v, float):
            print(f"{k}: {v:.4f}")
        else:
            print(f"{k}: {v}")

    print_confusion(
        "Class confusion",
        results["class_confusion"],
    )

    print_confusion(
        "Order confusion",
        results["order_confusion"],
    )

    print()


if __name__ == "__main__":
    main()
