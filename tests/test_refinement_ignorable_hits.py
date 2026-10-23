import json

import pandas as pd

from drayte.pipeline.runner import refinement_checkpoint_is_current
from drayte.refinement.refine_repeatmasker import (
    REFINEMENT_ALGORITHM_VERSION,
    RMHit,
    build_loci,
    refine_repeatmasker,
)


def make_hit(start, end, name, repeat_class):
    return RMHit(
        sw_score=100,
        perc_div=5.0,
        query="scaffold_1",
        q_start=start,
        q_end=end,
        strand="+",
        repeat_name=name,
        repeat_class=repeat_class,
        r_start="1",
        r_end=str(end - start + 1),
        r_left="(0)",
        rm_id=f"hit_{start}",
        length=end - start + 1,
        raw_index=start,
    )


def test_loose_mode_preserves_ignorable_hits_while_merging_te_fragments():
    hits = [
        make_hit(10, 20, "(AT)n", "Simple_repeat"),
        make_hit(100, 150, "TE_family_1", "DNA/TcMar"),
        make_hit(151, 160, "low_complexity", "Low_complexity"),
        make_hit(161, 210, "TE_family_1", "DNA/TcMar"),
        make_hit(220, 230, "5S", "rRNA"),
    ]

    loci = build_loci(
        hits=hits,
        species="Test",
        max_gap=150,
        mode="loose",
    )

    assert len(loci) == 4
    assert [locus.repeat_class for locus in loci] == [
        "Simple_repeat",
        "DNA/TcMar",
        "Low_complexity",
        "rRNA",
    ]
    te_locus = next(locus for locus in loci if locus.repeat_name == "TE_family_1")
    assert len(te_locus.hits) == 2
    assert (te_locus.start, te_locus.end) == (100, 210)


def test_old_refinement_manifest_is_not_reused(tmp_path):
    outputs = [tmp_path / name for name in ["refined.tsv", "filtered.tsv", "stats.tsv"]]
    for output in outputs:
        output.write_text("result\n")

    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"mode": "loose", "max_gap": 150}))

    assert not refinement_checkpoint_is_current(
        required_outputs=outputs,
        manifest_path=manifest,
        expected_parameters={"mode": "loose", "max_gap": 150},
    )

    manifest.write_text(
        json.dumps(
            {
                "refinement_algorithm_version": REFINEMENT_ALGORITHM_VERSION,
                "mode": "loose",
                "max_gap": 150,
            }
        )
    )

    assert refinement_checkpoint_is_current(
        required_outputs=outputs,
        manifest_path=manifest,
        expected_parameters={"mode": "loose", "max_gap": 150},
    )


def test_refinement_keeps_non_te_records_in_filtered_output(tmp_path):
    rmout = tmp_path / "test.fa.out"
    rmout.write_text(
        "   SW  perc div. perc del. perc ins. query sequence\n"
        "score  div. del. ins. sequence begin end (left) strand repeat class/family\n"
        "100 5.0 0.0 0.0 scaffold_1 10 20 (1000) + (AT)n Simple_repeat 1 11 (0) 1\n"
        "200 3.0 0.0 0.0 scaffold_1 100 150 (900) + TE_1 DNA/TcMar 1 51 (0) 2\n"
        "100 5.0 0.0 0.0 scaffold_1 200 210 (800) + low Low_complexity 1 11 (0) 3\n"
        "200 3.0 0.0 0.0 scaffold_1 300 350 (700) + TE_2 LINE/L1 1 51 (0) 4\n"
        "100 5.0 0.0 0.0 scaffold_1 400 410 (600) + 5S rRNA 1 11 (0) 5\n"
    )

    result = refine_repeatmasker(
        rmout=rmout,
        species="Test",
        outdir=tmp_path / "refinement",
        mode="loose",
        resolve_overlaps=True,
    )

    filtered = pd.read_csv(result["filtered_tsv"], sep="\t")
    assert result["raw_hits"] == 5
    assert result["preserved_ignorable_hits"] == 3
    assert set(filtered["class"]) == {
        "Simple_repeat",
        "Low_complexity",
        "rRNA",
        "DNA/TcMar",
        "LINE/L1",
    }
