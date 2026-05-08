from pathlib import Path

from drayte.classification.dfam import parse_nhmmer_tblout
from drayte.classification.dfam import (
    DfamHit,
    best_dfam_hits_by_family,
)

def test_parse_nhmmer_tblout(tmp_path: Path):
    tblout = tmp_path / "dfam.tblout"

    tblout.write_text(
        "# fake nhmmer tblout\n"
        "Helitron-4_Hbur#RC/Helitron - Helitron1Na_Mam DF000000006.4 "
        "702 906 389 196 409 176 2454 - 7.6e-10 34.5 4.9 @Heliconius\n"
        "weak_family - Arthur1A DF000000070.4 "
        "62 106 2288 2335 2271 2355 2418 + 0.83 7.9 11.3 @Heliconius\n"
    )

    hits = parse_nhmmer_tblout(tblout)

    assert len(hits) == 1
    assert hits[0].family_id == "Helitron-4_Hbur"
    assert hits[0].model_name == "Helitron1Na_Mam"
    assert hits[0].accession == "DF000000006.4"
    assert hits[0].evalue == 7.6e-10
    assert hits[0].score == 34.5

def test_best_dfam_hits_by_family():
    hits = [
        DfamHit(
            family_id="fam1",
            model_name="CR1_A",
            accession="DF1",
            evalue=1e-10,
            score=50,
            ali_start=1,
            ali_end=100,
        ),
        DfamHit(
            family_id="fam1",
            model_name="CR1_B",
            accession="DF2",
            evalue=1e-20,
            score=70,
            ali_start=1,
            ali_end=100,
        ),
        DfamHit(
            family_id="fam2",
            model_name="Gypsy_A",
            accession="DF3",
            evalue=1e-5,
            score=30,
            ali_start=1,
            ali_end=100,
        ),
    ]

    best = best_dfam_hits_by_family(hits)

    assert len(best) == 2
    assert best["fam1"].model_name == "CR1_B"
    assert best["fam2"].model_name == "Gypsy_A"


def test_split_fasta_round_robin(tmp_path: Path):
    from drayte.classification.dfam import split_fasta_round_robin

    fasta = tmp_path / "input.fa"
    fasta.write_text(
        ">a\nAAAA\n"
        ">b\nCCCC\n"
        ">c\nGGGG\n"
        ">d\nTTTT\n"
    )

    chunks = split_fasta_round_robin(
        fasta=fasta,
        outdir=tmp_path / "chunks",
        chunks=2,
    )

    assert len(chunks) == 2

    contents = [p.read_text() for p in chunks]

    assert ">a" in contents[0]
    assert ">c" in contents[0]
    assert ">b" in contents[1]
    assert ">d" in contents[1]


def test_merge_tblout_files(tmp_path: Path):
    from drayte.classification.dfam import merge_tblout_files

    a = tmp_path / "a.tblout"
    b = tmp_path / "b.tblout"
    out = tmp_path / "merged.tblout"

    a.write_text("# header a\nhit_a\n")
    b.write_text("# header b\nhit_b\n")

    merge_tblout_files([a, b], out)

    text = out.read_text()

    assert "hit_a" in text
    assert "hit_b" in text
