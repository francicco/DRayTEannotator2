from drayte.classification.ids import clean_family_id


def test_clean_family_id():

    assert (
        clean_family_id("Copia-2_Hm#LTR/Copia")
        == "Copia-2_Hm"
    )

    assert (
        clean_family_id("CR1-1_Hmel_A#LINE/CR1")
        == "CR1-1_Hmel_A"
    )

    assert (
        clean_family_id("SimpleFamily")
        == "SimpleFamily"
    )
