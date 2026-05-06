from drayte.classification.models import Family
from drayte.classification.classify import classify_family


def make_family(**kwargs):
    defaults = dict(
        family_id="test",
        consensus_len=1000,
        n_copies=1,
        homology_class="Unknown",
        homology_superfamily="Unknown",
        homology_score=0.0,
        domains=set(),
        ltr_present=False,
        tir_present=False,
        helitron_signal=False,
        tsd_present=False,
        polyA_present=False,
        orf_count=0,
        orf_max_len=0,
        boundary_consistency=0.0,
        fragmentation_score=1.0,
    )
    defaults.update(kwargs)
    return Family(**defaults)


def test_ltr_classification():
    f = make_family(
        homology_class="LTR",
        homology_superfamily="Gypsy",
        homology_score=0.8,
        domains={"RT_domain", "INT_domain"},
        ltr_present=True,
        tsd_present=True,
        boundary_consistency=0.9,
        fragmentation_score=0.1,
    )
    result = classify_family(f)
    assert result["class"] == "Class_I"
    assert result["order"] == "LTR"
    assert result["superfamily"] == "Gypsy"
    assert result["status"] == "OK"


def test_tir_classification():
    f = make_family(
        homology_class="DNA",
        homology_superfamily="hAT",
        homology_score=0.7,
        domains={"Transposase_domain"},
        tir_present=True,
        tsd_present=True,
        boundary_consistency=0.8,
        fragmentation_score=0.2,
    )
    result = classify_family(f)
    assert result["class"] == "Class_II"
    assert result["order"] == "TIR"
    assert result["superfamily"] == "hAT"


def test_line_classification():
    f = make_family(
        homology_class="LINE",
        homology_superfamily="L2",
        homology_score=0.6,
        domains={"RT_domain"},
        polyA_present=True,
        boundary_consistency=0.7,
        fragmentation_score=0.3,
    )
    result = classify_family(f)
    assert result["class"] == "Class_I"
    assert result["order"] == "LINE"
    assert result["superfamily"] == "L2"


def test_helitron_classification():
    f = make_family(
        homology_class="Helitron",
        homology_superfamily="Helitron",
        homology_score=0.85,
        helitron_signal=True,
        boundary_consistency=0.85,
        fragmentation_score=0.2,
    )
    result = classify_family(f)
    assert result["class"] == "Class_II"
    assert result["order"] == "Helitron"


def test_unknown_classification():
    f = make_family()
    result = classify_family(f)
    assert result["class"] == "Unknown"
    assert result["order"] == "Unknown"
    assert result["status"] == "unknown"
