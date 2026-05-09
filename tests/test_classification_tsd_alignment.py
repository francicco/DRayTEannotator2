from pathlib import Path

from drayte.classification.tsd_from_alignment import (
    detect_family_tsd_from_alignment,
)


def test_detect_family_tsd(tmp_path: Path):

    fa = tmp_path / "copies.fa"

    fa.write_text(
        ">copy1\n"
        "AAAATA" + "CCCCCCCC" + "TAGGGG\n"
        ">copy2\n"
        "GGGGTA" + "TTTTTTTT" + "TAGAAA\n"
        ">copy3\n"
        "CCCCTA" + "AAAAAAAA" + "TACCCC\n"
    )

    hit = detect_family_tsd_from_alignment(
        family_id="fam1",
        alignment_fasta=fa,
        flank_size=6,
        min_len=2,
        max_len=2,
    )

    assert hit.tsd_present is True
    assert hit.tsd_seq == "TA"
    assert hit.n_copies == 3
