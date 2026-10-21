import csv
import logging

from drayte.curation.organize import copy_extension_artifacts


def test_copy_extension_artifacts_matches_classified_directory_suffix(tmp_path):
    family_table = tmp_path / "family_table.tsv"
    with open(family_table, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["name", "top_orf_class"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "name": "Pmex-rnd-1_family-0",
                "top_orf_class": "DNA",
            }
        )

    extensionwork = tmp_path / "extensionwork"
    source_dir = extensionwork / "Pmex-rnd-1_family-0#Unknown"
    source_dir.mkdir(parents=True)
    (source_dir / "Pmex-rnd-1_family-0#Unknown_rep.fa").write_text(
        ">Pmex-rnd-1_family-0#Unknown\nACGT\n"
    )
    (source_dir / "Pmex-rnd-1_family-0#Unknown_MSA_extended.fa").write_text(
        ">Pmex-rnd-1_family-0#Unknown\nACGT\n"
    )
    (source_dir / "Pmex-rnd-1_family-0#Unknown.png").write_bytes(b"png")

    output_dirs = copy_extension_artifacts(
        family_table_tsv=family_table,
        extensionwork_dir=extensionwork,
        outdir=tmp_path / "te-aid",
        logger=logging.getLogger("test.curation-organize"),
    )

    dna_dir = output_dirs["DNA"]
    assert (dna_dir / "Pmex-rnd-1_family-0_rep.fa").exists()
    assert (dna_dir / "Pmex-rnd-1_family-0_MSA_extended.fa").exists()
    assert (dna_dir / "Pmex-rnd-1_family-0.png").exists()
