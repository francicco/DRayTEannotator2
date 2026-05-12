from pathlib import Path

from drayte.classification.ltr_detect import (
    detect_ltr_structure,
    detect_ltrs_from_fasta,
)


def make_ltr(length: int = 120) -> str:
    seed = (
        "ATGCCATTAGGCTAACCGTATCGGATCGTACGATTA"
        "CCGGAATTCGATCGGTTACGACTAGGCTAACGTTA"
    )
    return (seed * ((length // len(seed)) + 1))[:length]


def test_trim_like_terminal_direct_repeats():
    ltr = make_ltr(100)
    internal = "GATTACA" * 20
    seq = ltr + internal + ltr

    call = detect_ltr_structure("fam#Unknown", seq)

    assert call.family_id == "fam"
    assert call.ltr_present is True
    assert call.ltr_len == len(ltr)
    assert call.internal_len == len(internal)
    assert call.ltr_structural_type == "TRIM"


def test_lard_like_terminal_direct_repeats():
    ltr = make_ltr(120)
    internal = "GATTACA" * 300
    seq = ltr + internal + ltr

    call = detect_ltr_structure("fam", seq)

    assert call.ltr_present is True
    assert call.ltr_len == len(ltr)
    assert call.internal_len == len(internal)
    assert call.ltr_structural_type == "LARD_like"


def test_no_ltr_for_short_sequence():
    call = detect_ltr_structure("fam", "ACGT" * 20)

    assert call.ltr_present is False
    assert call.ltr_structural_type == "none"


def test_no_ltr_for_reverse_complement_repeat():
    left = make_ltr(100)
    table = str.maketrans("ACGT", "TGCA")
    right = left.translate(table)[::-1]
    seq = left + "GATTACA" * 20 + right

    call = detect_ltr_structure("fam", seq)

    assert call.ltr_present is False
    assert call.ltr_structural_type == "none"


def test_detect_ltrs_from_fasta_cleans_ids(tmp_path: Path):
    fasta = tmp_path / "lib.fa"
    ltr = make_ltr(100)
    fasta.write_text(">fam1#LTR/Unknown\n" + ltr + "GATTACA" * 20 + ltr + "\n")

    calls = detect_ltrs_from_fasta(fasta)

    assert len(calls) == 1
    assert calls[0].family_id == "fam1"
    assert calls[0].ltr_present is True
    assert calls[0].ltr_structural_type == "TRIM"
