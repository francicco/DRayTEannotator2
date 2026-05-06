from __future__ import annotations

import csv


def load_expected(path):
    expected = {}

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            expected[row["family_id"]] = row

    return expected


def evaluate_classifications(results, expected):

    metrics = {
        "total": 0,
        "correct_class": 0,
        "correct_order": 0,
        "correct_superfamily": 0,
    }

    for row in results:

        family_id = row["family_id"]

        if family_id not in expected:
            continue

        exp = expected[family_id]

        metrics["total"] += 1

        if row["class"] == exp["expected_class"]:
            metrics["correct_class"] += 1

        if row["order"] == exp["expected_order"]:
            metrics["correct_order"] += 1

        if row["superfamily"] == exp["expected_superfamily"]:
            metrics["correct_superfamily"] += 1

    return metrics
