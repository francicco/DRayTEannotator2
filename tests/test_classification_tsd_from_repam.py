from pathlib import Path

from drayte.classification.tsd_from_repam import (
    best_tsd_for_copy,
    infer_family_tsd_from_repam,
)


def test_best_tsd_for_copy():
    seq = "GGGTA" + "A" * 100 + "TAGGG"

    hit = best_tsd_for_copy(
        seq=seq,
        extended_left=5,
        extended_right=5,
        min_len=2,
        max_len=4,
    )

    assert hit is not None
    assert hit.tsd_len == 2
    assert hit.tsd_seq == "TA"


def test_infer_family_tsd_from_repam(tmp_path: Path):
    ranges = tmp_path / "repam-ranges.tsv"
    repseq = tmp_path / "repam-repseq.fa"

    ranges.write_text(
        "chr1\t1\t110\t+\tn=0,anchor_range=10-100,extended_left=5,extended_right=5\n"
        "chr1\t200\t309\t+\tn=1,anchor_range=10-100,extended_left=5,extended_right=5\n"
        "chr1\t400\t509\t+\tn=2,anchor_range=10-100,extended_left=5,extended_right=5\n"
        "chr1\t600\t709\t+\tn=3,anchor_range=10-100,extended_left=5,extended_right=5\n"
    )

    repseq.write_text(
        ">copy0\n"
        + "GGGTA" + "A" * 100 + "TAGGG" + "\n"
        ">copy1\n"
        + "CCCTA" + "C" * 100 + "TACCC" + "\n"
        ">copy2\n"
        + "TTTTA" + "G" * 100 + "TATTT" + "\n"
        ">copy3\n"
        + "GGGCC" + "T" * 100 + "AAGGG" + "\n"
    )

    call = infer_family_tsd_from_repam(
        repam_ranges=ranges,
        repam_repseq=repseq,
        family_id="fam1",
        min_len=2,
        max_len=4,
        min_copies=3,
        min_support=0.50,
    )

    assert call.family_id == "fam1"
    assert call.tsd_present is True
    assert call.tsd_len == 2
    assert call.tsd_seq == "TA"
    assert call.tsd_support == 0.75
