from pathlib import Path

from drayte.classification.pipeline import run_domain_annotation


def test_run_domain_annotation(tmp_path, monkeypatch):

    def fake_hmmscan(
        hmm_db,
        proteins_fasta,
        domtblout,
        hmmscan_bin="hmmscan",
        cpu=1,
    ):
        Path(domtblout).write_text(
            "# fake domtblout\n"
            "RVT_1 PF00078.1 200 fam1_orf1 - 300 1e-50 180.0 0.0 "
            "1 1 1e-50 1e-50 180.0 0.0 5 180 20 210 20 210 0.98 description\n"
        )
        return domtblout

    monkeypatch.setattr(
        "drayte.classification.pipeline.run_hmmscan",
        fake_hmmscan,
    )

    fasta = tmp_path / "families.fa"

    fasta.write_text(
        ">fam1\n"
        "ATG" + "AAA" * 200 + "TAA" + "\n"
    )

    hits = run_domain_annotation(
        consensus_fasta=fasta,
        hmm_db="fake.hmm",
        outdir=tmp_path / "domains",
        min_orf_nt=300,
        include_reverse_orfs=False,
    )

    assert len(hits) == 1
    assert hits[0].family_id == "fam1"
    assert hits[0].domain == "RVT_1"
