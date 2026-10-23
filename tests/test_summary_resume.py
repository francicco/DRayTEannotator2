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
            classification_result={
                "outdir": "/run/classification",
                "classified_library": "/run/classification/final-library.fa",
                "classifications_tsv": "/run/classification/classes.tsv",
                "evidence_tsv": "/run/classification/evidence.tsv",
                "nonauto_candidates_tsv": "/run/classification/nonauto.tsv",
                "n_families": 300,
                "n_classified": 237,
                "n_unknown": 63,
            },
            final_annotation_result={
                "outdir": "/run/final_annotation",
                "repeatmasker_out": "/run/final_annotation/final.fa.out",
                "repeatmasker_gff": "/run/final_annotation/final.fa.gff",
                "repeatmasker_tbl": "/run/final_annotation/final.fa.tbl",
                "masked_genome": "/run/final_annotation/final.fa.masked",
                "align": "/run/final_annotation/final.fa.align",
            },
            refinement_result={
                "outdir": "/run/refinement",
                "refined_tsv": "/run/refined.tsv",
                "refined_gff": "/run/refined.gff3",
                "refined_bed": "/run/refined.bed",
                "filtered_tsv": "/run/filtered.tsv",
                "filtered_gff": "/run/filtered.gff3",
                "filtered_bed": "/run/filtered.bed",
                "nested_gff": "/run/nested.gff3",
                "stats_tsv": "/run/stats.tsv",
            },
            summary_no_filter_result={
                "outdir": "/run/unfiltered_summary",
                "high_level_count": "/run/unfiltered_summary/high.txt",
                "family_level_count": "/run/unfiltered_summary/family.txt",
                "divergence_summary": "/run/unfiltered_summary/divergence.tsv",
                "defragmented_repeats": "/run/unfiltered_summary/repeats.tsv",
            },
            summary_filtered_result={
                "outdir": "/run/filtered_summary",
                "high_level_count": "/run/filtered_summary/high.txt",
                "family_level_count": "/run/filtered_summary/family.txt",
                "divergence_summary": "/run/filtered_summary/divergence.tsv",
                "defragmented_repeats": "/run/filtered_summary/repeats.tsv",
            },
            family_inspection_result={
                "manifest": "/run/families.tsv",
                "outdir": "/run/families",
                "families": 237,
            },
        )

    assert "FINAL RESULTS" in caplog.text
    assert "1. CLASSIFIED TE LIBRARY AND EVIDENCE" in caplog.text
    assert "Final classified library: /run/classification/final-library.fa" in caplog.text
    assert "2. FINAL REPEATMASKER GENOME ANNOTATION" in caplog.text
    assert "RepeatMasker GFF: /run/final_annotation/final.fa.gff" in caplog.text
    assert "3. REFINED AND FILTERED ANNOTATION" in caplog.text
    assert "Filtered non-overlapping loci (TSV): /run/filtered.tsv" in caplog.text
    assert "4. TE SUMMARY - UNFILTERED REFINED LOCI" in caplog.text
    assert "Plots and PDF summaries: /run/unfiltered_summary" in caplog.text
    assert "5. TE SUMMARY - FILTERED NON-OVERLAPPING LOCI" in caplog.text
    assert "Family-level count table: /run/filtered_summary/family.txt" in caplog.text
    assert "6. PER-FAMILY INSPECTION" in caplog.text
    assert "Family manifest: /run/families.tsv" in caplog.text
    assert "Families represented: 237" in caplog.text
