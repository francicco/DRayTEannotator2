from pathlib import Path

from drayte.classification.features import (
    consensus_lengths,
    build_families_from_evidence,
)
from drayte.classification.hmmer import DomainHit


def test_consensus_lengths(tmp_path: Path):
    fasta = tmp_path / "families.fa"
    fasta.write_text(">fam1\nATGC\n>fam2\nATGCAA\n")

    lengths = consensus_lengths(fasta)

    assert lengths["fam1"] == 4
    assert lengths["fam2"] == 6


def test_build_families_from_evidence(tmp_path: Path):
    fasta = tmp_path / "families.fa"
    fasta.write_text(
        ">fam1\n"
        "ATG" + "AAA" * 200 + "TAA" + "\n"
        ">fam2\n"
        "ATGAAATAA\n"
    )

    hits = [
        DomainHit(
            family_id="fam1",
            orf_id="fam1_orf1",
            domain="RVT_1",
            accession="PF00078",
            evalue=1e-50,
            score=180.0,
            hmm_start=1,
            hmm_end=100,
            ali_start=1,
            ali_end=100,
        )
    ]

    families = build_families_from_evidence(
        fasta,
        domain_hits=hits,
        min_orf_nt=300,
        include_reverse_orfs=False,
    )

    by_id = {f.family_id: f for f in families}

    assert by_id["fam1"].consensus_len > 300
    assert by_id["fam1"].orf_count >= 1
    assert by_id["fam1"].rt_present
    assert not by_id["fam2"].rt_present


def test_build_families_from_structure_evidence(tmp_path: Path):
    from drayte.classification.structure import StructureEvidence

    fasta = tmp_path / "families.fa"
    fasta.write_text(
        ">fam1\n"
        "ATG" + "AAA" * 200 + "TAA" + "\n"
    )

    families = build_families_from_evidence(
        fasta,
        domain_hits=[],
        structure_evidence=[
            StructureEvidence("fam1", "LTR", 0.9, "ltrharvest"),
            StructureEvidence("fam1", "TSD", 0.8, "ltrharvest"),
        ],
        min_orf_nt=300,
        include_reverse_orfs=False,
    )

    f = families[0]

    assert f.ltr_present
    assert f.tsd_present
    assert not f.tir_present


def test_build_families_from_dfam_evidence(tmp_path: Path):
    from drayte.classification.dfam import DfamHit

    fasta = tmp_path / "families.fa"
    fasta.write_text(
        ">fam1\n"
        "ATG" + "AAA" * 200 + "TAA" + "\n"
    )

    families = build_families_from_evidence(
        fasta,
        domain_hits=[],
        dfam_hits=[
            DfamHit(
                family_id="fam1",
                model_name="CR1-16_AMi",
                accession="DF000001265.3",
                evalue=1e-10,
                score=40.0,
                ali_start=1,
                ali_end=200,
            )
        ],
        min_orf_nt=300,
        include_reverse_orfs=False,
    )

    f = families[0]

    assert f.dfam_class == "Class_I"
    assert f.dfam_order == "LINE"
    assert f.dfam_superfamily == "CR1"
    assert f.dfam_model == "CR1-16_AMi"
    assert f.dfam_score == 40.0


def test_build_families_preserves_header_labels(tmp_path: Path):
    fasta = tmp_path / "families.fa"
    fasta.write_text(
        ">fam1#RC/Helitron @SpeciesX\n"
        "ATG" + "AAA" * 200 + "TAA" + "\n"
    )

    families = build_families_from_evidence(
        fasta,
        domain_hits=[],
        min_orf_nt=300,
        include_reverse_orfs=False,
    )

    f = families[0]

    assert f.family_id == "fam1"
    assert f.header_class == "RC"
    assert f.header_superfamily == "Helitron"
    assert f.original_header == "fam1#RC/Helitron"
