from pathlib import Path

from drayte.classification.sine_detect import (
    detect_sine_structure,
    detect_sines_from_fasta,
)
from drayte.classification.structure_merge import merge_structure_evidence
from drayte.classification.structure_detect import TirDetection


def test_polyA_tail_scores_sine_candidate():
    seq = "GCTAGCTAGC" * 20 + "A" * 18 + "TATATATATA"
    call = detect_sine_structure("fam#Unknown", seq)

    assert call.family_id == "fam"
    assert call.polyA_present is True
    assert call.polyA_tail_len >= 18
    assert call.sine_candidate is True
    assert call.sine_score >= 0.7


def test_sine_candidate_is_suppressed_by_tir():
    seq = "GCTAGCTAGC" * 20 + "A" * 20
    sine = detect_sine_structure("fam", seq)

    merged = merge_structure_evidence(
        tir_calls=[
            TirDetection(
                family_id="fam",
                tir_present=True,
                tir_len=22,
                tir_identity=0.95,
                confidence="HIGH",
            )
        ],
        tsd_calls=[],
        sine_calls=[sine],
        consensus_lengths={"fam": len(seq)},
    )

    assert merged[0].sine_candidate is False
    assert merged[0].sine_score == 0.0
    assert merged[0].structural_type.startswith("short_TIR")


def test_no_tir_sine_structural_type():
    seq = "GCTAGCTAGC" * 20 + "T" * 16
    sine = detect_sine_structure("fam", seq)

    merged = merge_structure_evidence(
        tir_calls=[],
        tsd_calls=[],
        sine_calls=[sine],
        consensus_lengths={"fam": len(seq)},
    )

    assert merged[0].sine_candidate is False
    assert merged[0].structural_type == "SINE_like"
    assert merged[0].polyA_present is True


def test_detect_sines_from_fasta_cleans_ids(tmp_path: Path):
    fasta = tmp_path / "lib.fa"
    fasta.write_text(">fam1#SINE/Unknown\n" + "ACGT" * 50 + "A" * 20 + "\n")

    calls = detect_sines_from_fasta(fasta)
    assert len(calls) == 1
    assert calls[0].family_id == "fam1"
    assert calls[0].sine_candidate is True
