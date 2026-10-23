from pathlib import Path

import drayte.classification.dfam as dfam_module
from drayte.classification.dfam import parse_nhmmer_tblout
from drayte.classification.dfam import (
    DfamHit,
    best_dfam_hits_by_family,
    nhmmer_tblout_complete,
    run_nhmmer_parallel_chunks,
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


def test_split_fasta_balanced_by_length(tmp_path: Path):
    from drayte.utils.fasta_split import split_fasta_balanced_by_length

    fasta = tmp_path / "input.fa"
    chunks_dir = tmp_path / "chunks"

    fasta.write_text(
        ">long\n" + "A" * 100 + "\n"
        ">medium\n" + "A" * 50 + "\n"
        ">short1\n" + "A" * 25 + "\n"
        ">short2\n" + "A" * 25 + "\n"
    )

    chunks = split_fasta_balanced_by_length(
        input_fasta=fasta,
        outdir=chunks_dir,
        n_chunks=2,
        prefix="chunk",
    )

    assert len(chunks) == 2
    assert all(p.exists() for p in chunks)

    sizes = []

    for path in chunks:
        total = 0
        current = []

        for line in path.read_text().splitlines():
            if line.startswith(">"):
                continue
            total += len(line.strip())

        sizes.append(total)

    assert sorted(sizes) == [100, 100]

    report = chunks_dir / "chunk.sizes.tsv"
    assert report.exists()

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


def test_nhmmer_tblout_complete_requires_ok_footer(tmp_path: Path):
    missing = tmp_path / "missing.tblout"
    empty = tmp_path / "empty.tblout"
    interrupted = tmp_path / "interrupted.tblout"
    complete = tmp_path / "complete.tblout"

    empty.write_text("")
    interrupted.write_text("# header\nhit\n")
    complete.write_text("# header\n# [ok]\n")

    assert not nhmmer_tblout_complete(missing)
    assert not nhmmer_tblout_complete(empty)
    assert not nhmmer_tblout_complete(interrupted)
    assert nhmmer_tblout_complete(complete)


def test_parallel_dfam_reuses_complete_legacy_chunks(
    tmp_path: Path,
    monkeypatch,
):
    consensus = tmp_path / "consensus.fa"
    dfam_db = tmp_path / "dfam.hmm"
    outdir = tmp_path / "dfam"
    chunks_dir = outdir / "chunks"
    tblout_dir = outdir / "tblouts"

    consensus.write_text(
        ">family_a\nAAAA\n"
        ">family_b\nCCCC\n"
    )
    dfam_db.write_text("HMM database\n")
    chunks = dfam_module.split_fasta_balanced_by_length(
        input_fasta=consensus,
        outdir=chunks_dir,
        n_chunks=2,
        prefix="chunk",
    )
    tblout_dir.mkdir(parents=True)
    (tblout_dir / f"{chunks[0].stem}.tblout").write_text(
        "# completed cached result\n# [ok]\n"
    )
    (tblout_dir / f"{chunks[1].stem}.tblout").write_text("")

    calls = []

    def fake_run_nhmmer(**kwargs):
        calls.append(Path(kwargs["consensus_fasta"]).name)
        result = Path(kwargs["tblout"])
        result.write_text("# newly completed result\n# [ok]\n")
        return result

    monkeypatch.setattr(dfam_module, "run_nhmmer", fake_run_nhmmer)

    merged = run_nhmmer_parallel_chunks(
        dfam_db=dfam_db,
        consensus_fasta=consensus,
        outdir=outdir,
        chunks=2,
        max_parallel=1,
    )

    assert calls == [chunks[1].name]
    assert merged.exists()
    assert merged.read_text().count("# [ok]") == 2
    assert (
        tblout_dir / f"{chunks[0].stem}.tblout.manifest.json"
    ).exists()
    assert (
        tblout_dir / f"{chunks[1].stem}.tblout.manifest.json"
    ).exists()


def test_parallel_dfam_reruns_only_chunk_whose_sequence_changed(
    tmp_path: Path,
    monkeypatch,
):
    consensus = tmp_path / "consensus.fa"
    dfam_db = tmp_path / "dfam.hmm"
    outdir = tmp_path / "dfam"
    dfam_db.write_text("HMM database\n")
    consensus.write_text(
        ">family_a\nAAAA\n"
        ">family_b\nCCCC\n"
    )

    calls = []

    def fake_run_nhmmer(**kwargs):
        chunk = Path(kwargs["consensus_fasta"])
        calls.append(chunk.read_text())
        result = Path(kwargs["tblout"])
        result.write_text("# completed result\n# [ok]\n")
        return result

    monkeypatch.setattr(dfam_module, "run_nhmmer", fake_run_nhmmer)

    kwargs = {
        "dfam_db": dfam_db,
        "consensus_fasta": consensus,
        "outdir": outdir,
        "chunks": 2,
        "max_parallel": 1,
    }
    run_nhmmer_parallel_chunks(**kwargs)
    assert len(calls) == 2

    run_nhmmer_parallel_chunks(**kwargs)
    assert len(calls) == 2

    consensus.write_text(
        ">family_a\nTTTT\n"
        ">family_b\nCCCC\n"
    )
    run_nhmmer_parallel_chunks(**kwargs)

    assert len(calls) == 3
    assert ">family_a\nTTTT\n" in calls[-1]
