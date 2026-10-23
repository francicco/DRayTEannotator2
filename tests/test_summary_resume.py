import logging

import pandas as pd

from drayte.pipeline.runner import log_final_results
from drayte.reporting.SummaryFilesGen import run_summary_files


def create_cached_summary(outdir, species):
    required_suffixes = [
        ".summaryPie.pdf",
        ".familyLevelCount.txt",
        ".highLevelCount.kable",
        ".familyLevelCount.kable",
        "_classification_landscape.pdf",
        "_split_class_landscape.pdf",
        "_superfamily_div_plot.pdf",
        "_divergence_summary_table.tsv",
        ".defragmented_repeats.tsv",
        ".highLevel_TE_composition.pdf",
    ]
    for suffix in required_suffixes:
        (outdir / f"{species}{suffix}").write_text("cached\n")

    pd.DataFrame(
        [
            {
                "TE Classification": "LINE",
                "Number of insertions": 42,
                "Coverage (bp)": 12345,
                "% Genome Coverage": 1.25,
            }
        ]
    ).to_csv(outdir / f"{species}.highLevelCount.txt", sep="\t", index=False)


def test_cached_summary_is_skipped_but_printed_and_returns_full_contract(
    tmp_path,
    capsys,
    caplog,
):
    species = "Mmim"
    outdir = tmp_path / "summary_Filtered"
    outdir.mkdir()
    create_cached_summary(outdir, species)

    with caplog.at_level(logging.INFO):
        result = run_summary_files(
            refined_tsv=tmp_path / "cached-input.tsv",
            genome=tmp_path / "genome.fa",
            species=species,
            outdir=outdir,
        )

    printed = capsys.readouterr().out
    assert "Summary files already exist; skipping" in caplog.text
    assert f"{species} TE summary (high-level)" in printed
    assert f"Results directory: {outdir}" in printed
    assert "LINE" in printed
    assert "12345" in printed
    assert result["high_level_count"] == str(outdir / f"{species}.highLevelCount.txt")
    assert result["family_level_count"] == str(outdir / f"{species}.familyLevelCount.txt")


def test_final_results_are_logged_after_resumed_stages(caplog):
    logger = logging.getLogger("drayte-test")

    with caplog.at_level(logging.INFO, logger="drayte-test"):
        log_final_results(
            logger=logger,
            classification_result={"classified_library": "/run/final-library.fa"},
            final_annotation_result={"repeatmasker_out": "/run/final.fa.out"},
            refinement_result={
                "refined_tsv": "/run/refined.tsv",
                "filtered_tsv": "/run/filtered.tsv",
                "stats_tsv": "/run/stats.tsv",
            },
            summary_no_filter_result={"high_level_count": "/run/unfiltered.txt"},
            summary_filtered_result={"high_level_count": "/run/filtered.txt"},
            family_inspection_result={
                "manifest": "/run/families.tsv",
                "outdir": "/run/families",
                "families": 237,
            },
        )

    assert "FINAL RESULTS" in caplog.text
    assert "Classified TE library: /run/final-library.fa" in caplog.text
    assert "Filtered high-level TE summary: /run/filtered.txt" in caplog.text
    assert "Family inspection: /run/families.tsv (237 families)" in caplog.text
