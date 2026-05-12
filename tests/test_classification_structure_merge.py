from drayte.classification.complexity_detect import ComplexityDetection
from drayte.classification.sine_detect import SineDetection
from drayte.classification.structure_detect import TirDetection
from drayte.classification.structure_merge import merge_structure_evidence
from drayte.classification.tsd_from_repam import RepamTSDCall


def make_tir(
    family_id: str,
    tir_len: int = 25,
    tir_identity: float = 0.95,
    confidence: str = "HIGH",
    **kwargs,
) -> TirDetection:
    return TirDetection(
        family_id=family_id,
        tir_present=True,
        tir_len=tir_len,
        tir_identity=tir_identity,
        confidence=confidence,
        **kwargs,
    )


def make_tsd(
    family_id: str,
    tsd_len: int = 2,
    tsd_seq: str = "TA",
    tsd_support: float = 0.70,
) -> RepamTSDCall:
    return RepamTSDCall(
        family_id=family_id,
        tsd_present=True,
        tsd_len=tsd_len,
        tsd_seq=tsd_seq,
        tsd_support=tsd_support,
        n_copies_checked=10,
        n_copies_with_tsd=7,
        n_distinct_tsd_hits=1,
    )


def make_simple_repeat(family_id: str, consensus_len: int = 300) -> ComplexityDetection:
    return ComplexityDetection(
        family_id=family_id,
        consensus_len=consensus_len,
        seq_entropy=0.4,
        dominant_1mer_fraction=0.85,
        dominant_2mer_fraction=0.92,
        dominant_3mer_fraction=0.2,
        dominant_4mer_fraction=0.2,
        max_homopolymer_run=20,
        simple_repeat_candidate=True,
        low_complexity_candidate=False,
        complexity_class="simple_repeat",
        complexity_confidence="HIGH",
    )


def make_low_complexity(family_id: str, consensus_len: int = 300) -> ComplexityDetection:
    return ComplexityDetection(
        family_id=family_id,
        consensus_len=consensus_len,
        seq_entropy=1.1,
        dominant_1mer_fraction=0.72,
        dominant_2mer_fraction=0.40,
        dominant_3mer_fraction=0.35,
        dominant_4mer_fraction=0.35,
        max_homopolymer_run=8,
        simple_repeat_candidate=False,
        low_complexity_candidate=True,
        complexity_class="low_complexity",
        complexity_confidence="MEDIUM",
    )


def make_sine(family_id: str, consensus_len: int = 250) -> SineDetection:
    return SineDetection(
        family_id=family_id,
        consensus_len=consensus_len,
        polyA_present=True,
        polyA_tail_len=20,
        polyA_fraction=0.8,
        tail_type="polyA",
        poliii_a_box=True,
        poliii_b_box=True,
        poliii_score=1.0,
        sine_head_type="tRNA_like",
        sine_head_score=0.7,
        sine_score=1.0,
        sine_candidate=True,
        confidence="HIGH",
    )


def test_merge_structure_evidence():
    merged = merge_structure_evidence(
        tir_calls=[
            make_tir(
                family_id="fam1",
                tir_len=22,
                tir_identity=0.94,
            )
        ],
        tsd_calls=[
            make_tsd(
                family_id="fam1",
                tsd_support=0.55,
            )
        ],
        consensus_lengths={"fam1": 450},
    )

    assert len(merged) == 1

    x = merged[0]
    assert x.structure_class == "TIR_TSD"
    assert x.structural_type == "MITE"
    assert x.tir_grade == "STRONG"
    assert x.superfamily_hint == "Tc1-Mariner_like"
    assert x.superfamily_confidence == "HIGH"
    assert x.mite_candidate is True
    assert x.tir_tsd_element is True
    assert x.confidence in {"HIGH", "MEDIUM"}


def test_tir_superfamily_hints():
    examples = [
        (2, "TA", "Tc1-Mariner_like", "HIGH"),
        (4, "TTAA", "piggyBac_like", "HIGH"),
        (3, "TAA", "PIF-Harbinger_like", "MEDIUM"),
        (2, "CA", "2bp_TSD_TIR_ambiguous", "LOW"),
        (8, "ACGTACGT", "hAT_P_Merlin_ambiguous", "LOW"),
        (10, "ACGTACGTAA", "Mutator_like", "MEDIUM"),
    ]

    for tsd_len, tsd_seq, hint, conf in examples:
        merged = merge_structure_evidence(
            tir_calls=[make_tir("fam")],
            tsd_calls=[make_tsd("fam", tsd_len=tsd_len, tsd_seq=tsd_seq, tsd_support=0.55)],
            consensus_lengths={"fam": 450},
        )

        assert merged[0].superfamily_hint == hint
        assert merged[0].superfamily_confidence == conf


def test_cacta_hint_requires_terminal_motif_for_2bp_non_ta_tsd():
    merged = merge_structure_evidence(
        tir_calls=[
            make_tir(
                family_id="fam_cacta",
                left_tir_seq="CACTAGGGTTTAAACCCGGGTTTAA",
                right_tir_seq="TTAAACCCGGGTTTAAACCCTAGTG",
            )
        ],
        tsd_calls=[make_tsd("fam_cacta", tsd_len=2, tsd_seq="CA", tsd_support=0.55)],
        consensus_lengths={"fam_cacta": 450},
    )

    assert merged[0].tir_5p_motif.startswith("CACTA")
    assert merged[0].tir_3p_motif.startswith("CACTA")
    assert merged[0].superfamily_hint == "CACTA_like"
    assert merged[0].superfamily_confidence == "MEDIUM"


def test_terminal_motifs_are_written_to_structure_evidence():
    merged = merge_structure_evidence(
        tir_calls=[
            make_tir(
                family_id="fam_motif",
                tir_len=20,
                tir_identity=1.0,
                left_tir_seq="TACGAAACCCGGG",
                right_tir_seq="CCCGGGTTTCGTA",
            )
        ],
        tsd_calls=[],
        consensus_lengths={"fam_motif": 500},
    )

    assert merged[0].tir_5p_motif == "TACGAAAC"
    assert merged[0].tir_3p_motif == "TACGAAAC"
    assert merged[0].tir_terminal_motif == "TACGAAAC"


def test_long_tir_tsd_element_is_not_called_mite_like():
    merged = merge_structure_evidence(
        tir_calls=[make_tir("long_tir_tsd", tir_len=30, tir_identity=0.96)],
        tsd_calls=[make_tsd("long_tir_tsd")],
        consensus_lengths={"long_tir_tsd": 1500},
    )

    x = merged[0]
    assert x.structural_type == "TIR_TSD_element"
    assert x.tir_tsd_element is True
    assert x.mite_candidate is False


def test_short_weak_tir_tsd_is_not_strict_mite():
    merged = merge_structure_evidence(
        tir_calls=[
            make_tir(
                family_id="short_weak",
                tir_len=10,
                tir_identity=0.80,
                confidence="LOW",
            )
        ],
        tsd_calls=[make_tsd("short_weak", tsd_len=2, tsd_seq="CA", tsd_support=0.45)],
        consensus_lengths={"short_weak": 350},
    )

    x = merged[0]
    assert x.structural_type == "short_TIR_TSD_element"
    assert x.tir_tsd_element is True
    assert x.mite_candidate is False


def test_simple_repeat_fallback():
    merged = merge_structure_evidence(
        tir_calls=[],
        tsd_calls=[],
        complexity_calls=[make_simple_repeat("sr")],
        consensus_lengths={"sr": 300},
    )

    x = merged[0]
    assert x.structural_type == "simple_repeat"
    assert x.simple_repeat_candidate is True
    assert x.low_complexity_candidate is False
    assert x.complexity_class == "simple_repeat"


def test_low_complexity_fallback():
    merged = merge_structure_evidence(
        tir_calls=[],
        tsd_calls=[],
        complexity_calls=[make_low_complexity("lc")],
        consensus_lengths={"lc": 300},
    )

    x = merged[0]
    assert x.structural_type == "low_complexity"
    assert x.simple_repeat_candidate is False
    assert x.low_complexity_candidate is True
    assert x.complexity_class == "low_complexity"


def test_tir_overrides_low_complexity():
    merged = merge_structure_evidence(
        tir_calls=[make_tir("tir_lc")],
        tsd_calls=[],
        complexity_calls=[make_low_complexity("tir_lc", consensus_len=400)],
        consensus_lengths={"tir_lc": 400},
    )

    x = merged[0]
    assert x.structural_type == "short_TIR_element"
    assert x.low_complexity_candidate is True


def test_tir_tsd_overrides_simple_repeat():
    merged = merge_structure_evidence(
        tir_calls=[make_tir("tir_sr")],
        tsd_calls=[make_tsd("tir_sr")],
        complexity_calls=[make_simple_repeat("tir_sr", consensus_len=500)],
        consensus_lengths={"tir_sr": 500},
    )

    x = merged[0]
    assert x.structural_type == "MITE"
    assert x.simple_repeat_candidate is True


def test_sine_overrides_low_complexity():
    merged = merge_structure_evidence(
        tir_calls=[],
        tsd_calls=[],
        sine_calls=[make_sine("sine_lc")],
        complexity_calls=[make_low_complexity("sine_lc", consensus_len=250)],
        consensus_lengths={"sine_lc": 250},
    )

    x = merged[0]
    assert x.structural_type == "SINE"
    assert x.low_complexity_candidate is True
