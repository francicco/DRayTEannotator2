from drayte.classification.parsers.repeatmasker import (
    parse_repeatmasker_header,
    classification_from_repeatmasker,
)


def test_parse_repeatmasker_header():

    parsed = parse_repeatmasker_header(
        "Gypsy1#LTR/Gypsy"
    )

    assert parsed.family_id == "Gypsy1"
    assert parsed.rm_class == "LTR"
    assert parsed.rm_superfamily == "Gypsy"


def test_classification_from_repeatmasker():

    parsed = parse_repeatmasker_header(
        "Helitron1#RC/Helitron"
    )

    result = classification_from_repeatmasker(parsed)

    assert result["expected_class"] == "Class_II"
    assert result["expected_order"] == "Helitron"
    assert result["expected_superfamily"] == "Helitron"
