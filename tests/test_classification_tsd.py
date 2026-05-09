from drayte.classification.tsd import (
    best_tsd_from_boundaries,
    detect_tsd_for_copy,
    summarize_tsd_hits,
)


def test_detect_simple_tsd():
    hit = best_tsd_from_boundaries(
        left_flank="AAAAATC",
        right_flank="ATCGGGG",
    )

    assert hit is not None

    seq, identity = hit

    assert seq == "ATC"
    assert identity == 1.0


def test_no_tsd():
    hit = best_tsd_from_boundaries(
        left_flank="AAAAAAA",
        right_flank="GGGGGGG",
    )

    assert hit is None


def test_detect_tsd_for_copy():
    hit = detect_tsd_for_copy(
        family_id="fam1",
        left_flank="CCCCCTTAA",
        right_flank="TTAAGGGGG",
    )

    assert hit.tsd_present is True
    assert hit.tsd_seq == "TTAA"
    assert hit.tsd_len == 4


def test_summarize_tsd_hits():
    hits = [
        detect_tsd_for_copy("fam1", "AAAATA", "TAGGGG"),
        detect_tsd_for_copy("fam1", "CCCCTA", "TAGGGG"),
        detect_tsd_for_copy("fam1", "GGGGTA", "TAGGGG"),
        detect_tsd_for_copy("fam1", "AAAAAA", "CCCCCC"),
    ]

    summary = summarize_tsd_hits("fam1", hits)

    assert summary.tsd_present is True
    assert summary.tsd_seq == "TA"
    assert summary.n_copies == 4
