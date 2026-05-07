from drayte.classification.dfammap import infer_te_from_dfam_model


def test_dfam_helitron():
    assert infer_te_from_dfam_model("Helitron1Na_Mam") == (
        "Class_II",
        "Helitron",
        "Helitron",
    )


def test_dfam_cr1():
    assert infer_te_from_dfam_model("CR1-16_AMi") == (
        "Class_I",
        "LINE",
        "CR1",
    )


def test_dfam_rte():
    assert infer_te_from_dfam_model("MamRTE1") == (
        "Class_I",
        "LINE",
        "RTE",
    )


def test_dfam_mariner():
    assert infer_te_from_dfam_model("Mariner3_CE") == (
        "Class_II",
        "TIR",
        "TcMar-Mariner",
    )
