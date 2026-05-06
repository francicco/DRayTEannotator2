from drayte.classification.benchmark import (
    evaluate_classifications,
)


def test_evaluate_classifications():

    results = [
        {
            "family_id": "fam1",
            "class": "Class_I",
            "order": "LTR",
            "superfamily": "Gypsy",
        },
        {
            "family_id": "fam2",
            "class": "Class_II",
            "order": "TIR",
            "superfamily": "hAT",
        },
    ]

    expected = {
        "fam1": {
            "expected_class": "Class_I",
            "expected_order": "LTR",
            "expected_superfamily": "Gypsy",
        },
        "fam2": {
            "expected_class": "Class_II",
            "expected_order": "LINE",
            "expected_superfamily": "hAT",
        },
    }

    metrics = evaluate_classifications(
        results,
        expected,
    )

    assert metrics["total"] == 2
    assert metrics["correct_class"] == 2
    assert metrics["correct_order"] == 1
    assert metrics["correct_superfamily"] == 2
