from argparse import Namespace
from pathlib import Path

import pandas as pd

from compare_earlgrey_drayann import run


def write_summary(path: Path, total_percent: float, dna_percent: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    genome_size = 1000
    rows = [
        {
            "TE Classification": "DNA",
            "Coverage (bp)": int(dna_percent * 10),
            "Copy Number": 5,
            "% Genome Coverage": dna_percent,
            "Genome Size": genome_size,
            "TE Family Count": 2,
        },
        {
            "TE Classification": "DNA-nested",
            "Coverage (bp)": 10,
            "Copy Number": 1,
            "% Genome Coverage": 1.0,
            "Genome Size": genome_size,
            "TE Family Count": 1,
        },
        {
            "TE Classification": "Unclassified",
            "Coverage (bp)": int((total_percent - dna_percent - 1) * 10),
            "Copy Number": 3,
            "% Genome Coverage": total_percent - dna_percent - 1,
            "Genome Size": genome_size,
            "TE Family Count": 1,
        },
        {
            "TE Classification": "Total Interspersed Repeat",
            "Coverage (bp)": int(total_percent * 10),
            "Copy Number": 9,
            "% Genome Coverage": total_percent,
            "Genome Size": genome_size,
            "TE Family Count": pd.NA,
        },
        {
            "TE Classification": "Non-Repeat",
            "Coverage (bp)": int((100 - total_percent) * 10),
            "Copy Number": "notApplicable",
            "% Genome Coverage": 100 - total_percent,
            "Genome Size": genome_size,
            "TE Family Count": pd.NA,
        },
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_compare_discovers_pairs_and_collapses_nested(tmp_path: Path):
    root = tmp_path / "annotations"
    outdir = tmp_path / "comparison"
    earl = root / "Sp1_Earlgrey" / "Sp1_summaryFiles" / "Sp1.highLevelCount.txt"
    dray = (
        root
        / "Sp1_DRayTEannotator2"
        / "Sp1.Defragmenting_RModAnnotation_Filtered"
        / "Sp1.highLevelCount.txt"
    )
    write_summary(earl, total_percent=20.0, dna_percent=5.0)
    write_summary(dray, total_percent=25.0, dna_percent=7.0)

    result = run(
        Namespace(
            root=root,
            outdir=outdir,
            drayann_summary="filtered",
            species=None,
            no_plots=False,
        )
    )

    assert result == 0
    comparison = pd.read_csv(outdir / "class_comparison.tsv", sep="\t")
    dna = comparison[comparison["TE Classification"] == "DNA"].iloc[0]
    assert dna["% Genome Coverage EarlGrey"] == 6.0
    assert dna["% Genome Coverage DRayAnn"] == 8.0
    assert dna["% Genome Coverage difference (DRayAnn - EarlGrey)"] == 2.0

    totals = pd.read_csv(outdir / "species_total_comparison.tsv", sep="\t")
    assert totals.iloc[0]["Total repeat coverage (%) EarlGrey"] == 20.0
    assert totals.iloc[0]["Total repeat coverage (%) DRayAnn"] == 25.0
    assert (outdir / "comparison_report.md").exists()
    assert (outdir / "total_repeat_coverage_filtered.png").exists()
    assert (outdir / "total_repeat_coverage_scatter_filtered.pdf").exists()
    assert (outdir / "class_coverage_difference_heatmap_filtered.png").exists()


def test_compare_reports_missing_partner(tmp_path: Path):
    root = tmp_path / "annotations"
    outdir = tmp_path / "comparison"
    earl = root / "Sp1_EarlGrey" / "Sp1_summaryFiles" / "Sp1.highLevelCount.txt"
    write_summary(earl, total_percent=20.0, dna_percent=5.0)

    result = run(
        Namespace(
            root=root,
            outdir=outdir,
            drayann_summary="filtered",
            species=None,
            no_plots=True,
        )
    )

    assert result == 0
    inventory = pd.read_csv(outdir / "input_inventory.tsv", sep="\t")
    missing = inventory[inventory["Method"] == "DRayAnn filtered"].iloc[0]
    assert missing["Status"] == "missing"
    comparison = pd.read_csv(outdir / "species_total_comparison.tsv", sep="\t")
    assert comparison.empty
