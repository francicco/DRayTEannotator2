#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from drayte.utils.paths import ensure_dir, stage_dir


CLASS_ORDER = [
    "DNA", "DNA-nested",
    "Rolling Circle", "Rolling Circle-nested",
    "Penelope", "Penelope-nested",
    "LINE", "LINE-nested",
    "SINE", "SINE-nested",
    "LTR", "LTR-nested",
    "Other (Simple Repeat, Microsatellite, RNA)",
    "Other (Simple Repeat, Microsatellite, RNA)-nested",
    "Unclassified", "Unclassified-nested",
]


def fasta_size(fasta: Path) -> int:
    total = 0
    with fasta.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def major_class(repeat_class: str) -> str:
    x = str(repeat_class)

    if not x or x == "nan":
        return "Unclassified"

    if "Penelope" in x or x.startswith("PLE"):
        return "Penelope"
    if x.startswith("DNA"):
        return "DNA"
    if x.startswith("LINE"):
        return "LINE"
    if x.startswith("SINE"):
        return "SINE"
    if x.startswith("LTR"):
        return "LTR"
    if x.startswith("RC") or "Helitron" in x or "Rolling" in x:
        return "Rolling Circle"
    if any(
        k in x
        for k in [
            "Simple_repeat",
            "Low_complexity",
            "Satellite",
            "RNA",
            "rRNA",
            "tRNA",
            "snRNA",
            "srpRNA",
        ]
    ):
        return "Other (Simple Repeat, Microsatellite, RNA)"
    if "Unknown" in x or "Unclassified" in x:
        return "Unclassified"

    return "Unclassified"


def read_repeatmasker_out(path: Path) -> pd.DataFrame:
    rows = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue
            if line.startswith(("SW", "score", "perc", "---")):
                continue

            parts = line.split()

            if len(parts) < 14:
                continue

            try:
                sw_score = float(parts[0])
                perc_div = float(parts[1])
                perc_del = float(parts[2])
                perc_ins = float(parts[3])
                query = parts[4]
                start = int(parts[5])
                end = int(parts[6])
                strand = parts[8]
                family = parts[9]
                class_family = parts[10]
            except ValueError:
                continue

            if "/" in class_family:
                repeat_class, superfamily = class_family.split("/", 1)
            else:
                repeat_class = class_family
                superfamily = class_family

            cls = major_class(class_family)

            rows.append(
                {
                    "query": query,
                    "start": min(start, end),
                    "end": max(start, end),
                    "length_bp": abs(end - start) + 1,
                    "strand": strand,
                    "family": family,
                    "repeat_class_raw": class_family,
                    "repeat_class": repeat_class,
                    "class": cls,
                    "classification": cls,
                    "superfamily": superfamily or "Unknown",
                    "perc_div": perc_div,
                    "perc_del": perc_del,
                    "perc_ins": perc_ins,
                    "sw_score": sw_score,
                    "repeatmasker_nested_flag": parts[-1] == "*",
                }
            )

    df = pd.DataFrame(rows)

    if df.empty:
        raise RuntimeError(f"No RepeatMasker records parsed from {path}")

    return df


def defragment_hits(df: pd.DataFrame, max_gap: int = 100) -> pd.DataFrame:
    merged_rows = []

    sort_cols = ["query", "family", "repeat_class_raw", "strand", "start", "end"]
    df = df.sort_values(sort_cols).reset_index(drop=True)

    group_cols = ["query", "family", "repeat_class_raw", "strand"]

    for _, group in df.groupby(group_cols, dropna=False, sort=False):
        group = group.sort_values(["start", "end"])

        current = None
        div_values = []
        score_values = []

        for row in group.to_dict("records"):
            if current is None:
                current = row.copy()
                div_values = [row["perc_div"]]
                score_values = [row["sw_score"]]
                current["fragments_merged"] = 1
                continue

            gap = row["start"] - current["end"] - 1

            if gap <= max_gap:
                old_len = current["length_bp"]
                new_len = row["length_bp"]

                current["end"] = max(current["end"], row["end"])
                current["start"] = min(current["start"], row["start"])
                current["length_bp"] = current["end"] - current["start"] + 1
                current["fragments_merged"] += 1

                div_values.append(row["perc_div"])
                score_values.append(row["sw_score"])

                total_len = old_len + new_len
                current["perc_div"] = (
                    (current["perc_div"] * old_len) + (row["perc_div"] * new_len)
                ) / total_len
                current["perc_del"] = (
                    (current["perc_del"] * old_len) + (row["perc_del"] * new_len)
                ) / total_len
                current["perc_ins"] = (
                    (current["perc_ins"] * old_len) + (row["perc_ins"] * new_len)
                ) / total_len
                current["sw_score"] = sum(score_values)
                current["repeatmasker_nested_flag"] = (
                    current["repeatmasker_nested_flag"] or row["repeatmasker_nested_flag"]
                )
            else:
                merged_rows.append(current)
                current = row.copy()
                div_values = [row["perc_div"]]
                score_values = [row["sw_score"]]
                current["fragments_merged"] = 1

        if current is not None:
            merged_rows.append(current)

    out = pd.DataFrame(merged_rows)
    out["length_bp"] = out["end"] - out["start"] + 1
    return out.reset_index(drop=True)


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def detect_nested_hits(
    df: pd.DataFrame,
    min_overlap_fraction: float = 0.80,
    min_outer_extra_bp: int = 100,
) -> pd.DataFrame:
    df = df.copy()
    df["nested"] = False

    nested_flags = []

    for query, group in df.groupby("query", sort=False):
        group = group.sort_values(["start", "end"]).reset_index()

        for _, hit in group.iterrows():
            hit_len = hit["length_bp"]
            is_nested = False

            candidates = group[
                (group["start"] <= hit["start"])
                & (group["end"] >= hit["end"])
                & (group["index"] != hit["index"])
            ]

            for _, outer in candidates.iterrows():
                if outer["family"] == hit["family"]:
                    continue

                outer_len = outer["length_bp"]

                if outer_len < hit_len + min_outer_extra_bp:
                    continue

                overlap = interval_overlap(
                    int(hit["start"]),
                    int(hit["end"]),
                    int(outer["start"]),
                    int(outer["end"]),
                )

                if overlap / hit_len >= min_overlap_fraction:
                    is_nested = True
                    break

            nested_flags.append((hit["index"], is_nested))

    for idx, flag in nested_flags:
        df.loc[idx, "nested"] = flag

    df["classification"] = df.apply(
        lambda r: f"{r['class']}-nested" if r["nested"] else r["class"],
        axis=1,
    )

    return df


def clean_repeatmasker_hits(
    df: pd.DataFrame,
    max_merge_gap: int = 100,
    min_nested_overlap_fraction: float = 0.80,
) -> pd.DataFrame:
    df = defragment_hits(df, max_gap=max_merge_gap)
    df = detect_nested_hits(df, min_overlap_fraction=min_nested_overlap_fraction)
    return df


def write_markdown_table(df: pd.DataFrame, path: Path) -> None:
    path.write_text(df.to_markdown(index=False) + "\n")


def make_tables(
    df: pd.DataFrame,
    genome_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    family = (
        df.groupby(
            ["classification", "family", "repeat_class_raw", "class", "superfamily"],
            dropna=False,
        )
        .agg(
            **{
                "Coverage (bp)": ("length_bp", "sum"),
                "Copy Number": ("length_bp", "size"),
                "Mean Divergence": ("perc_div", "mean"),
                "Median Divergence": ("perc_div", "median"),
            }
        )
        .reset_index()
    )

    family["% Genome Coverage"] = 100 * family["Coverage (bp)"] / genome_size
    family["Genome Size"] = genome_size
    family = family.sort_values("Coverage (bp)", ascending=False)

    family_out = family.rename(
        columns={
            "classification": "TE Classification",
            "family": "Family",
            "repeat_class_raw": "RepeatMasker Classification",
            "superfamily": "Superfamily",
        }
    )[
        [
            "TE Classification",
            "Family",
            "RepeatMasker Classification",
            "Superfamily",
            "Coverage (bp)",
            "Copy Number",
            "% Genome Coverage",
            "Genome Size",
            "Mean Divergence",
            "Median Divergence",
        ]
    ]

    high = (
        family.groupby("classification", dropna=False)
        .agg(
            **{
                "Coverage (bp)": ("Coverage (bp)", "sum"),
                "Copy Number": ("Copy Number", "sum"),
                "% Genome Coverage": ("% Genome Coverage", "sum"),
                "Genome Size": ("Genome Size", "first"),
                "TE Family Count": ("family", "nunique"),
            }
        )
        .reset_index()
        .rename(columns={"classification": "TE Classification"})
    )

    missing = []
    present = set(high["TE Classification"])

    for cls in CLASS_ORDER:
        if cls not in present:
            missing.append(
                {
                    "TE Classification": cls,
                    "Coverage (bp)": 0,
                    "Copy Number": 0,
                    "% Genome Coverage": 0.0,
                    "Genome Size": genome_size,
                    "TE Family Count": pd.NA,
                }
            )

    if missing:
        high = pd.concat([high, pd.DataFrame(missing)], ignore_index=True)

    high["TE Classification"] = pd.Categorical(
        high["TE Classification"],
        categories=CLASS_ORDER,
        ordered=True,
    )
    high = high.sort_values("TE Classification")
    high["TE Classification"] = high["TE Classification"].astype(str)

    total_bp = int(high["Coverage (bp)"].sum())
    total_copy = int(high["Copy Number"].fillna(0).sum())
    non_repeat_bp = genome_size - total_bp

    high_final = pd.concat(
        [
            high,
            pd.DataFrame(
                [
                    {
                        "TE Classification": "Total Interspersed Repeat",
                        "Coverage (bp)": total_bp,
                        "Copy Number": str(total_copy),
                        "% Genome Coverage": 100 * total_bp / genome_size,
                        "Genome Size": genome_size,
                        "TE Family Count": pd.NA,
                    },
                    {
                        "TE Classification": "Non-Repeat",
                        "Coverage (bp)": non_repeat_bp,
                        "Copy Number": "notApplicable",
                        "% Genome Coverage": 100 * non_repeat_bp / genome_size,
                        "Genome Size": genome_size,
                        "TE Family Count": pd.NA,
                    },
                ]
            ),
        ],
        ignore_index=True,
    )

    div = (
        df.assign(divergence_bin=df["perc_div"].astype(int))
        .groupby(["classification", "class", "superfamily", "divergence_bin"], dropna=False)
        .agg(
            coverage_bp=("length_bp", "sum"),
            copy_number=("length_bp", "size"),
        )
        .reset_index()
    )

    div["genome_percent"] = 100 * div["coverage_bp"] / genome_size

    return high_final, family_out, div


def save_pie(high: pd.DataFrame, species: str, outpath: Path) -> None:
    dat = high[
        ~high["TE Classification"].isin(["Total Interspersed Repeat", "Non-Repeat"])
    ].copy()
    dat = dat[dat["Coverage (bp)"] > 0]

    plt.figure(figsize=(8, 8))
    plt.pie(dat["Coverage (bp)"], labels=dat["TE Classification"], autopct="%1.1f%%")
    plt.title(f"{species} repeat composition")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def save_class_landscape(div: pd.DataFrame, species: str, outpath: Path) -> None:
    pivot = div.pivot_table(
        index="divergence_bin",
        columns="classification",
        values="genome_percent",
        aggfunc="sum",
        fill_value=0,
    ).sort_index()

    ax = pivot.plot(kind="bar", stacked=True, figsize=(12, 6), width=1.0)
    ax.set_xlabel("RepeatMasker divergence (%)")
    ax.set_ylabel("% genome coverage")
    ax.set_title(f"{species} repeat landscape")
    ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def save_split_landscape(div: pd.DataFrame, species: str, outpath: Path) -> None:
    classes = [c for c in CLASS_ORDER if c in set(div["classification"])]
    n = len(classes)

    if n == 0:
        return

    ncols = 2
    nrows = max(1, (n + 1) // 2)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(12, max(5, nrows * 2.3)),
    )
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for ax, cls in zip(axes, classes):
        sub = div[div["classification"] == cls]
        agg = sub.groupby("divergence_bin")["genome_percent"].sum()
        ax.bar(agg.index, agg.values, width=1.0)
        ax.set_title(cls, fontsize=9)
        ax.set_xlabel("Divergence (%)")
        ax.set_ylabel("% genome")

    for ax in axes[len(classes):]:
        ax.axis("off")

    fig.suptitle(f"{species} repeat landscape by classification")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def save_superfamily_plot(
    div: pd.DataFrame,
    species: str,
    outpath: Path,
    top_n: int = 20,
) -> None:
    top = (
        div.groupby("superfamily")["coverage_bp"]
        .sum()
        .sort_values(ascending=False)
        .head(top_n)
        .index
    )

    if len(top) == 0:
        return

    sub = div[div["superfamily"].isin(top)]

    ncols = 4
    nrows = max(1, (len(top) + ncols - 1) // ncols)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(14, max(6, nrows * 2.2)),
    )
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for ax, sf in zip(axes, top):
        sfdat = sub[sub["superfamily"] == sf]
        agg = sfdat.groupby("divergence_bin")["genome_percent"].sum()
        ax.bar(agg.index, agg.values, width=1.0)
        ax.set_title(str(sf), fontsize=8)
        ax.set_xlabel("Div. (%)")
        ax.set_ylabel("% genome")

    for ax in axes[len(top):]:
        ax.axis("off")

    fig.suptitle(f"{species} top superfamily divergence profiles")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def summary_files_exist(outdir: Path, species: str) -> bool:
    required = [
        outdir / f"{species}.summaryPie.pdf",
        outdir / f"{species}.highLevelCount.txt",
        outdir / f"{species}.familyLevelCount.txt",
        outdir / f"{species}.highLevelCount.kable",
        outdir / f"{species}.familyLevelCount.kable",
        outdir / f"{species}_classification_landscape.pdf",
        outdir / f"{species}_split_class_landscape.pdf",
        outdir / f"{species}_superfamily_div_plot.pdf",
        outdir / f"{species}_divergence_summary_table.tsv",
        outdir / f"{species}.defragmented_repeats.tsv",
    ]

    return all(f.exists() and f.stat().st_size > 0 for f in required)


def run_summary(
    rmout: Path,
    genome_size: int,
    species: str,
    outdir: Path,
    max_merge_gap: int = 100,
    min_nested_overlap_fraction: float = 0.80,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    raw = read_repeatmasker_out(rmout)
    clean = clean_repeatmasker_hits(
        raw,
        max_merge_gap=max_merge_gap,
        min_nested_overlap_fraction=min_nested_overlap_fraction,
    )

    clean.to_csv(outdir / f"{species}.defragmented_repeats.tsv", sep="\t", index=False)

    high, family, div = make_tables(clean, genome_size)

    high.to_csv(outdir / f"{species}.highLevelCount.txt", sep="\t", index=False)
    family.to_csv(outdir / f"{species}.familyLevelCount.txt", sep="\t", index=False)
    div.to_csv(outdir / f"{species}_divergence_summary_table.tsv", sep="\t", index=False)

    write_markdown_table(high, outdir / f"{species}.highLevelCount.kable")
    write_markdown_table(family, outdir / f"{species}.familyLevelCount.kable")

    save_pie(high, species, outdir / f"{species}.summaryPie.pdf")
    save_class_landscape(div, species, outdir / f"{species}_classification_landscape.pdf")
    save_split_landscape(div, species, outdir / f"{species}_split_class_landscape.pdf")
    save_superfamily_plot(div, species, outdir / f"{species}_superfamily_div_plot.pdf")


def run_summary_files(config, final_annotation_result: dict, logger) -> dict:
    species = config.species
    outdir = ensure_dir(stage_dir(config.outdir_path, "summaryFiles"))

    logger.info("=" * 80)
    logger.info("STAGE: summaryFiles")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    if summary_files_exist(outdir, species):
        logger.info("Summary files already exist; skipping")
        return {
            "stage": "summaryFiles",
            "outdir": str(outdir),
        }

    rmout = Path(final_annotation_result["repeatmasker_out"])
    if not rmout.exists() or rmout.stat().st_size == 0:
        raise FileNotFoundError(f"RepeatMasker .out file not found or empty: {rmout}")

    genome_fa = config.outdir_path / "discovery" / "assemblies_dir" / f"{species}.fa"
    if not genome_fa.exists() or genome_fa.stat().st_size == 0:
        raise FileNotFoundError(f"Genome FASTA not found or empty: {genome_fa}")

    genome_size = fasta_size(genome_fa)

    max_merge_gap = int(config.extra.get("summary_max_merge_gap", 100))
    min_nested_overlap_fraction = float(
        config.extra.get("summary_min_nested_overlap_fraction", 0.80)
    )

    logger.info("Genome size inferred from %s: %s bp", genome_fa, genome_size)
    logger.info("Generating summary files from %s", rmout)
    logger.info("Defragmentation max gap: %s bp", max_merge_gap)
    logger.info("Nested overlap fraction: %.2f", min_nested_overlap_fraction)

    run_summary(
        rmout=rmout,
        genome_size=genome_size,
        species=species,
        outdir=outdir,
        max_merge_gap=max_merge_gap,
        min_nested_overlap_fraction=min_nested_overlap_fraction,
    )

    logger.info("Summary files generated")

    return {
        "stage": "summaryFiles",
        "outdir": str(outdir),
        "high_level_count": str(outdir / f"{species}.highLevelCount.txt"),
        "family_level_count": str(outdir / f"{species}.familyLevelCount.txt"),
        "divergence_summary": str(outdir / f"{species}_divergence_summary_table.tsv"),
        "defragmented_repeats": str(outdir / f"{species}.defragmented_repeats.tsv"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate EarlGrey-like repeat summary files.")
    parser.add_argument("--rmout", required=True, type=Path)
    parser.add_argument("--genome-size", required=True, type=int)
    parser.add_argument("--species", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--max-merge-gap", type=int, default=100)
    parser.add_argument("--min-nested-overlap-fraction", type=float, default=0.80)

    args = parser.parse_args()

    run_summary(
        rmout=args.rmout,
        genome_size=args.genome_size,
        species=args.species,
        outdir=args.outdir,
        max_merge_gap=args.max_merge_gap,
        min_nested_overlap_fraction=args.min_nested_overlap_fraction,
    )


if __name__ == "__main__":
    main()
