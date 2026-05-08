from pathlib import Path

from drayte.classification.mmseqs import (
    best_mmseqs_hits_by_family,
    infer_rescue_evidence,
    parse_mmseqs_tsv,
)

def test_parse_mmseqs_tsv_and_rescue(tmp_path: Path):
    tsv = tmp_path / "hits.tsv"

    tsv.write_text(
        "fam1#Unknown\tfam1#Unknown\t1.000\t1000\t0\t0\t1\t1000\t1\t1000\t0\t2000\n"
        "fam1#Unknown\tfam2#DNA/hAT\t0.910\t900\t10\t0\t1\t900\t1\t900\t1e-50\t1800\n"
        "fam2#DNA/hAT\tfam1#Unknown\t0.910\t900\t10\t0\t1\t900\t1\t900\t1e-50\t1800\n"
        "fam1#Unknown\tfam3#Unknown\t0.950\t950\t5\t0\t1\t950\t1\t950\t1e-60\t1900\n"
    )

    hits = parse_mmseqs_tsv(tsv)

    assert len(hits) == 3

    best = best_mmseqs_hits_by_family(hits)
    assert best["fam1"].target_id == "fam3"

    rescue = infer_rescue_evidence(hits)
    assert "fam1" in rescue
    assert rescue["fam1"].rescue_order == "TIR"
    assert rescue["fam1"].rescue_superfamily == "Unknown"
