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


def test_infer_family_tsd_uses_modal_length_when_exact_sequence_varies(tmp_path: Path):
    ranges = tmp_path / "repam-ranges.tsv"
    repseq = tmp_path / "repam-repseq.fa"

    ranges.write_text(
        "chr1\t1\t110\t+\tn=0,anchor_range=10-100,extended_left=5,extended_right=5\n"
        "chr1\t200\t309\t+\tn=1,anchor_range=10-100,extended_left=5,extended_right=5\n"
        "chr1\t400\t509\t+\tn=2,anchor_range=10-100,extended_left=5,extended_right=5\n"
        "chr1\t600\t709\t+\tn=3,anchor_range=10-100,extended_left=5,extended_right=5\n"
        "chr1\t800\t909\t+\tn=4,anchor_range=10-100,extended_left=5,extended_right=5\n"
    )

    # Four copies have a 2-bp TSD, but the exact sequence differs.
    # This should still support a length-based TSD call for superfamily hints.
    repseq.write_text(
        ">copy0\n" + "GGGTA" + "A" * 100 + "TAGGG" + "\n"
        ">copy1\n" + "CCCTG" + "C" * 100 + "TGCCC" + "\n"
        ">copy2\n" + "TTTCA" + "G" * 100 + "CATTT" + "\n"
        ">copy3\n" + "AAACG" + "T" * 100 + "CGAAA" + "\n"
        ">copy4\n" + "GGGAA" + "T" * 100 + "CCGGG" + "\n"
    )

    call = infer_family_tsd_from_repam(
        repam_ranges=ranges,
        repam_repseq=repseq,
        family_id="fam_len",
        min_len=2,
        max_len=4,
        min_copies=3,
        min_support=0.50,
    )

    assert call.tsd_present is True
    assert call.tsd_len == 2
    assert call.tsd_support == 0.8
    assert call.tsd_len_support == 0.8
    assert call.tsd_exact_support == 0.2
    assert call.tsd_confidence in {"MEDIUM", "HIGH"}


def test_best_tsd_for_copy_tolerates_boundary_slop():
    # True TSD is TA. The provided left boundary is shifted one base to the
    # right, so a strict boundary-only search would miss it.
    seq = "GGGTA" + "A" * 100 + "TAGGG"

    hit = best_tsd_for_copy(
        seq=seq,
        extended_left=6,
        extended_right=5,
        min_len=2,
        max_len=4,
        boundary_slop=2,
    )

    assert hit is not None
    assert hit.tsd_len == 2
    assert hit.tsd_seq == "TA"
    assert hit.left_offset == -1
