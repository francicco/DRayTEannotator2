from drayte.classification.structure_detect import TirDetection
from drayte.classification.structure_merge import (
    merge_structure_evidence,
)


from drayte.classification.tsd_from_repam import (
    RepamTSDCall,
)


def test_merge_structure_evidence():
    tir_calls = [
        TirDetection(
            family_id="fam1",
            tir_present=True,
            tir_len=22,
            tir_identity=0.94,
            confidence="HIGH",
        )
    ]

    tsd_calls = [
        RepamTSDCall(
            family_id="fam1",
            tsd_present=True,
            tsd_len=2,
            tsd_seq="TA",
            tsd_support=0.55,
            n_copies_checked=10,
            n_copies_with_tsd=7,
            n_distinct_tsd_hits=1,
        )
    ]

    merged = merge_structure_evidence(
        tir_calls=tir_calls,
        tsd_calls=tsd_calls,
        consensus_lengths={
            "fam1": 450
        },
    )

    assert len(merged) == 1

    x = merged[0]

    assert x.structure_class == "TIR_TSD"
    assert x.mite_candidate is True
    assert x.confidence in {"HIGH", "MEDIUM"}
