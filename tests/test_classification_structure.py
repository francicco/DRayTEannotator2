from drayte.classification.structure import (
    StructureEvidence,
    normalize_structure_type,
    summarize_structure_evidence,
)


def test_normalize_structure_type():
    assert normalize_structure_type("LTR_retrotransposon") == "LTR"
    assert normalize_structure_type("TIR_structure") == "TIR"
    assert normalize_structure_type("Rolling_circle") == "HELITRON"
    assert normalize_structure_type("target_site_duplication") == "TSD"
    assert normalize_structure_type("polyA_tail") == "POLYA"
    assert normalize_structure_type("unknown") is None


def test_summarize_structure_evidence():
    evidence = [
        StructureEvidence("fam1", "LTR", 0.9, "ltrharvest"),
        StructureEvidence("fam1", "TSD", 0.8, "ltrharvest"),
        StructureEvidence("fam2", "HELITRON", 0.7, "heliano"),
        StructureEvidence("fam3", "unknown", 1.0, "test"),
    ]

    summary = summarize_structure_evidence(evidence)

    assert summary["fam1"]["ltr_present"]
    assert summary["fam1"]["tsd_present"]
    assert not summary["fam1"]["tir_present"]
    assert summary["fam1"]["best_structure_score"] == 0.9
    assert summary["fam1"]["structure_sources"] == ["ltrharvest"]

    assert summary["fam2"]["helitron_signal"]
    assert "fam3" not in summary
