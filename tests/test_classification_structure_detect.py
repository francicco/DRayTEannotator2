from pathlib import Path

from drayte.classification.structure_detect import (
    best_terminal_inverted_repeat,
    detect_tirs_from_fasta,
    write_tir_structure_tsv,
)


def test_best_terminal_inverted_repeat_detects_simple_tir():
    left_tir = "ATGCGTACGTTAGCTA"
    right_tir = "TAGCTAACGTACGCAT"  # reverse complement of left_tir
    seq = left_tir + ("A" * 100) + right_tir

    det = best_terminal_inverted_repeat(
        seq,
        window=50,
        min_len=15,
        min_identity=0.80,
    )

    assert det is not None
    assert det.tir_present
    assert det.tir_len >= 15
    assert det.tir_identity >= 0.80


def test_detect_tirs_from_fasta_and_write_tsv(tmp_path: Path):
    fasta = tmp_path / "families.fa"
    out = tmp_path / "structure.tsv"

    left_tir = "ATGCGTACGTTAGCTA"
    right_tir = "TAGCTAACGTACGCAT"

    fasta.write_text(
        ">fam1#DNA/hAT\n"
        + left_tir + ("A" * 100) + right_tir + "\n"
        ">fam2\n"
        + "A" * 200 + "\n"
    )

    detections = detect_tirs_from_fasta(
        fasta,
        window=50,
        min_len=15,
        min_identity=0.80,
    )

    by_id = {d.family_id: d for d in detections}

    assert by_id["fam1"].tir_present
    assert not by_id["fam2"].tir_present

    write_tir_structure_tsv(detections, out)

    text = out.read_text()
    assert "fam1\tTIR" in text
    assert "fam2" not in text
