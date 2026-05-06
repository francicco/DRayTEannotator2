import subprocess
import sys
from pathlib import Path


def test_classification_cli_tsv_mode(tmp_path: Path):
    families = tmp_path / "families.tsv"
    out = tmp_path / "classifications.tsv"

    families.write_text(
        "family_id\tconsensus_len\tn_copies\thomology_class\thomology_superfamily\thomology_score\t"
        "rt_present\tintegrase_present\ttransposase_present\tltr_present\ttir_present\thelitron_signal\t"
        "tsd_present\tpolyA_present\torf_count\torf_max_len\tboundary_consistency\tfragmentation_score\n"
        "fam1\t6000\t20\tLTR\tGypsy\t0.8\tTrue\tTrue\tFalse\tTrue\tFalse\tFalse\tTrue\tFalse\t2\t3000\t0.9\t0.1\n"
    )

    subprocess.run(
        [
            sys.executable,
            "-m",
            "drayte.classification.run",
            "--input",
            str(families),
            "--output",
            str(out),
        ],
        check=True,
    )

    text = out.read_text()
    assert "fam1\tClass_I\tLTR\tGypsy\tOK\tHIGH" in text


def test_classification_cli_fasta_mode(tmp_path: Path):
    fasta = tmp_path / "families.fa"
    out = tmp_path / "classifications.tsv"

    fasta.write_text(
        ">fam1\n"
        "ATG" + "AAA" * 200 + "TAA" + "\n"
    )

    subprocess.run(
        [
            sys.executable,
            "-m",
            "drayte.classification.run",
            "--fasta",
            str(fasta),
            "--output",
            str(out),
            "--min-orf-nt",
            "300",
            "--forward-only",
        ],
        check=True,
    )

    text = out.read_text()
    assert text.startswith("family_id\tclass\torder")
    assert "fam1\tUnknown\tUnknown\tUnknown" in text
