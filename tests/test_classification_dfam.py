from pathlib import Path

from drayte.classification.dfam import parse_nhmmer_tblout


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
