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


def test_evidence_from_structure_candidates(tmp_path):
    from drayte.structure.models import StructureCandidate
    from drayte.classification.structure import evidence_from_structure_candidates

    candidates = [
        StructureCandidate(
            candidate_id="fam_ltr",
            module="LTR",
            contig="ctg1",
            start=1,
            end=1000,
            strand="+",
            length=1000,
            fasta_path=tmp_path / "fam_ltr.fa",
        )
    ]

    evidence = evidence_from_structure_candidates(candidates)

    assert len(evidence) == 1
    assert evidence[0].family_id == "fam_ltr"
    assert evidence[0].evidence_type == "LTR"
    assert evidence[0].source == "structure:LTR"


def test_structure_evidence_tsv_io(tmp_path):
    from drayte.classification.structure import (
        StructureEvidence,
        write_structure_evidence_tsv,
        load_structure_evidence_tsv,
    )

    out = tmp_path / "structure.tsv"

    evidence = [
        StructureEvidence(
            family_id="fam1",
            evidence_type="LTR",
            score=0.9,
            source="ltrharvest",
        ),
        StructureEvidence(
            family_id="fam2",
            evidence_type="HELITRON",
            score=0.8,
            source="heliano",
        ),
    ]

    write_structure_evidence_tsv(
        evidence,
        out,
    )

    loaded = load_structure_evidence_tsv(out)

    assert len(loaded) == 2
    assert loaded[0].family_id == "fam1"
    assert loaded[0].evidence_type == "LTR"
    assert loaded[1].evidence_type == "HELITRON"
