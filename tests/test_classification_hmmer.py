from pathlib import Path

from drayte.classification.hmmer import (
    family_from_orf_id,
    parse_domtblout,
    summarize_domains_by_family,
)


def test_family_from_orf_id():
    assert family_from_orf_id("fam1_orf2") == "fam1"
    assert family_from_orf_id("fam1") == "fam1"


def test_parse_domtblout_and_summarize(tmp_path: Path):
    domtbl = tmp_path / "test.domtblout"

    domtbl.write_text(
        "# fake hmmscan domtblout\n"
        "RT_domain PF00078.1 200 fam1_orf1 - 300 1e-50 180.0 0.0 "
        "1 1 1e-50 1e-50 180.0 0.0 5 180 20 210 20 210 0.98 description\n"
        "INT_domain PF00665.1 250 fam1_orf1 - 300 1e-20 90.0 0.0 "
        "1 1 1e-20 1e-20 90.0 0.0 10 220 30 250 30 250 0.95 description\n"
        "Weak_domain PF00000.1 100 fam2_orf1 - 200 1e-1 5.0 0.0 "
        "1 1 1e-1 1e-1 5.0 0.0 1 50 1 50 1 50 0.50 description\n"
    )

    hits = parse_domtblout(domtbl, max_evalue=1e-5, min_score=20.0)

    assert len(hits) == 2
    assert hits[0].family_id == "fam1"
    assert hits[0].domain == "RT_domain"

    summary = summarize_domains_by_family(hits)

    assert "fam1" in summary
    assert summary["fam1"]["domain_hits"] == 2
    assert summary["fam1"]["best_domain_score"] == 180.0
    assert summary["fam1"]["domains"] == ["INT_domain", "RT_domain"]
    assert "fam2" not in summary
