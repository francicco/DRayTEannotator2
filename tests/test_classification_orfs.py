from pathlib import Path

from drayte.classification.orfs import (
    extract_orf_calls_by_family,
    summarize_orfs,
    write_translated_orfs,
)


def test_orf_summary_and_translation(tmp_path: Path):
    fasta = tmp_path / "families.fa"
    outpep = tmp_path / "orfs.aa.fa"

    fasta.write_text(
        ">fam1\n"
        "ATG" + "AAA" * 200 + "TAA" + "\n"
        ">fam2\n"
        "ATGAAATAA\n"
    )

    calls = extract_orf_calls_by_family(
        fasta,
        minsize_nt=300,
        include_reverse=False,
    )

    summary = summarize_orfs(calls)

    assert summary["fam1"]["orf_count"] >= 1
    assert summary["fam1"]["orf_max_len"] >= 300
    assert summary["fam2"]["orf_count"] == 0
    assert summary["fam2"]["orf_max_len"] == 0

    n = write_translated_orfs(
        fasta,
        outpep,
        minsize_nt=300,
        include_reverse=False,
    )

    assert n >= 1
    assert outpep.exists()
    assert outpep.stat().st_size > 0
