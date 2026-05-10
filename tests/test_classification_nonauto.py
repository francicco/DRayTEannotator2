from drayte.classification.models import Family
from drayte.classification.nonauto import (
    score_sine_candidate,
    score_mite_candidate,
    score_trim_candidate,
)


def test_sine_candidate_short_noncoding_polya():
    f = Family(
        family_id="sine1",
        consensus_len=300,
        n_copies=20,
        polyA_present=True,
        orf_count=0,
        orf_max_len=0,
        domains=set(),
    )

    call = score_sine_candidate(f)

    assert call is not None
    assert call.candidate_type == "SINE_candidate"


def test_mite_candidate_short_tir_noncoding():
    f = Family(
        family_id="mite1",
        consensus_len=450,
        n_copies=30,
        tir_present=True,
        tsd_present=True,
        orf_count=0,
        orf_max_len=0,
        domains=set(),
    )

    call = score_mite_candidate(f)

    assert call is not None
    assert call.candidate_type == "MITE_candidate"


def test_trim_candidate_short_ltr_noncoding():
    f = Family(
        family_id="trim1",
        consensus_len=900,
        n_copies=15,
        ltr_present=True,
        orf_count=0,
        orf_max_len=0,
        domains=set(),
    )

    call = score_trim_candidate(f)

    assert call is not None
    assert call.candidate_type == "TRIM_candidate"


def test_autonomous_line_is_not_sine_candidate():
    f = Family(
        family_id="line1",
        consensus_len=6000,
        n_copies=10,
        polyA_present=True,
        orf_count=1,
        orf_max_len=3000,
        domains={"RT"},
    )

    assert score_sine_candidate(f) is None


def test_sine_excludes_tir():
    f = Family(
        family_id="sine_like_but_tir",
        consensus_len=300,
        n_copies=20,
        polyA_present=True,
        tir_present=True,
        orf_count=0,
        orf_max_len=0,
        domains=set(),
    )

    assert score_sine_candidate(f) is None


def test_mite_tsd_superfamily_hint():
    f = Family(
        family_id="mite_ta",
        consensus_len=450,
        n_copies=30,
        tir_present=True,
        tsd_present=True,
        tsd_len=2,
        tsd_seq="TA",
        tsd_support=0.6,
        orf_count=0,
        orf_max_len=0,
        domains=set(),
    )

    call = score_mite_candidate(f)

    assert call is not None
    assert "tsd_superfamily_hint=TcMar-Mariner" in call.evidence
