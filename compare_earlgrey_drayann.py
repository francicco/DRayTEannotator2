#!/usr/bin/env python3
"""Compare final Earl Grey and DRayAnn annotation summaries across species."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


NUMERIC_COLUMNS = (
    "Coverage (bp)",
    "Copy Number",
    "% Genome Coverage",
    "Genome Size",
    "TE Family Count",
)
AGGREGATE_CLASSES = {"Total Interspersed Repeat", "Non-Repeat"}
CLASS_ORDER = (
    "DNA",
    "Rolling Circle",
    "Penelope",
    "LINE",
    "SINE",
    "LTR",
    "Other (Simple Repeat, Microsatellite, RNA)",
    "Unclassified",
    "Total Interspersed Repeat",
    "Non-Repeat",
)
CLASS_COLORS = {
    "DNA": "#4C78A8",
    "Rolling Circle": "#72B7B2",
    "Penelope": "#54A24B",
    "LINE": "#F58518",
    "SINE": "#ECA82C",
    "LTR": "#E45756",
    "Other (Simple Repeat, Microsatellite, RNA)": "#B279A2",
    "Unclassified": "#79706E",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare Earl Grey summaryFiles with DRayAnn final annotation "
            "summaries across species."
        )
    )
    parser.add_argument(
        "--root",
        required=True,
        type=Path,
        help="Directory containing SPECIES_Earlgrey and SPECIES_DRayTEannotator2 folders.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: ROOT/EarlGrey_vs_DRayAnn).",
    )
    parser.add_argument(
        "--drayann-summary",
        choices=("filtered", "unfiltered", "both"),
        default="filtered",
        help=(
            "DRayAnn result to compare with Earl Grey. The filtered final "
            "annotation is the default."
        ),
    )
    parser.add_argument(
        "--species",
        nargs="*",
        default=None,
        help="Optional species prefixes to include (for example Mmim Mmex).",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Write tables and report without PNG/PDF plots.",
    )
    return parser.parse_args(argv)


def method_directories(root: Path) -> dict[str, dict[str, Path]]:
    discovered: dict[str, dict[str, Path]] = {}
    for path in sorted(root.iterdir()):
        if not path.is_dir():
            continue
        lowered = path.name.casefold()
        if lowered.endswith("_earlgrey"):
            species = path.name[: -len("_Earlgrey")]
            discovered.setdefault(species.casefold(), {})["EarlGrey"] = path
            discovered[species.casefold()]["species"] = species
        elif lowered.endswith("_drayteannotator2"):
            species = path.name[: -len("_DRayTEannotator2")]
            discovered.setdefault(species.casefold(), {})["DRayAnn"] = path
            discovered[species.casefold()]["species"] = species
    return discovered


def first_existing(candidates: list[Path]) -> Path | None:
    for candidate in candidates:
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate
    return None


def find_earlgrey_summary(directory: Path, species: str) -> Path | None:
    direct = first_existing(
        [
            directory / f"{species}_summaryFiles" / f"{species}.highLevelCount.txt",
            directory / "summaryFiles" / f"{species}.highLevelCount.txt",
        ]
    )
    if direct is not None:
        return direct
    matches = sorted(directory.glob("*_summaryFiles/*.highLevelCount.txt"))
    return first_existing(matches)


def find_drayann_summary(
    directory: Path,
    species: str,
    summary_kind: str,
) -> Path | None:
    suffix = (
        "Defragmenting_RModAnnotation_Filtered"
        if summary_kind == "filtered"
        else "Defragmenting_RModAnnotation_noFiltered"
    )
    direct = directory / f"{species}.{suffix}" / f"{species}.highLevelCount.txt"
    if direct.exists() and direct.stat().st_size > 0:
        return direct

    wanted_parent = suffix.casefold()
    matches = [
        path
        for path in directory.rglob("*.highLevelCount.txt")
        if path.parent.name.casefold().endswith(wanted_parent)
    ]
    return first_existing(sorted(matches))


def normalize_classification(value: object) -> str:
    text = re.sub(r"\s+", " ", str(value).strip())
    nested = text.casefold().endswith("-nested")
    base = re.sub(r"-nested$", "", text, flags=re.IGNORECASE).strip()
    folded = base.casefold()

    if folded in {"rolling-circle", "rolling circle", "rc/helitron"}:
        base = "Rolling Circle"
    elif folded in {"unknown", "unclassified"}:
        base = "Unclassified"
    elif folded.startswith("other") and any(
        token in folded for token in ("simple", "microsatellite", "rna")
    ):
        base = "Other (Simple Repeat, Microsatellite, RNA)"
    elif folded == "total interspersed repeat":
        base = "Total Interspersed Repeat"
    elif folded in {"non-repeat", "non repeat"}:
        base = "Non-Repeat"

    return f"{base}-nested" if nested else base


def read_high_level(path: Path, species: str, method: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str)
    table.columns = [str(column).strip() for column in table.columns]
    required = {"TE Classification", "Coverage (bp)", "% Genome Coverage"}
    missing = sorted(required.difference(table.columns))
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(missing)}")

    for column in NUMERIC_COLUMNS:
        if column not in table.columns:
            table[column] = pd.NA
        table[column] = pd.to_numeric(table[column], errors="coerce")

    table["TE Classification"] = table["TE Classification"].map(
        normalize_classification
    )
    table["Nested"] = table["TE Classification"].str.endswith("-nested")
    table["Base Classification"] = table["TE Classification"].str.replace(
        r"-nested$", "", regex=True
    )
    table.insert(0, "Method", method)
    table.insert(0, "Species", species)
    table["Source file"] = str(path)
    return table


def collapse_nested_rows(raw: pd.DataFrame) -> pd.DataFrame:
    class_rows = raw[~raw["TE Classification"].isin(AGGREGATE_CLASSES)].copy()
    collapsed = (
        class_rows.groupby(
            ["Species", "Method", "Base Classification"],
            dropna=False,
            as_index=False,
        )
        .agg(
            {
                "Coverage (bp)": "sum",
                "Copy Number": "sum",
                "% Genome Coverage": "sum",
                "Genome Size": "max",
                # Parent and nested rows may refer to the same families.
                "TE Family Count": "max",
                "Source file": "first",
            }
        )
        .rename(columns={"Base Classification": "TE Classification"})
    )
    collapsed["Nested"] = False
    collapsed["Base Classification"] = collapsed["TE Classification"]

    aggregates = raw[raw["TE Classification"].isin(AGGREGATE_CLASSES)].copy()
    return pd.concat([collapsed, aggregates], ignore_index=True, sort=False)


def safe_ratio(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    return numerator.div(denominator.where(denominator != 0))


def pair_with_earlgrey(
    collapsed: pd.DataFrame,
    drayann_method: str,
) -> pd.DataFrame:
    keys = ["Species", "TE Classification"]
    metrics = list(NUMERIC_COLUMNS)
    earl = collapsed[collapsed["Method"] == "EarlGrey"][keys + metrics]
    dray = collapsed[collapsed["Method"] == drayann_method][keys + metrics]
    paired = earl.merge(
        dray,
        on=keys,
        how="inner",
        suffixes=(" EarlGrey", " DRayAnn"),
    )
    paired.insert(1, "DRayAnn summary", drayann_method.removeprefix("DRayAnn "))

    for metric in metrics:
        earl_column = f"{metric} EarlGrey"
        dray_column = f"{metric} DRayAnn"
        paired[f"{metric} difference (DRayAnn - EarlGrey)"] = (
            paired[dray_column] - paired[earl_column]
        )
        paired[f"{metric} ratio (DRayAnn / EarlGrey)"] = safe_ratio(
            paired[dray_column], paired[earl_column]
        )
    return paired


def species_summary(collapsed: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for (species, method), group in collapsed.groupby(["Species", "Method"]):
        indexed = group.set_index("TE Classification")

        def value(classification: str, column: str) -> float:
            if classification not in indexed.index:
                return 0.0
            result = indexed.loc[classification, column]
            if isinstance(result, pd.Series):
                result = result.iloc[0]
            return 0.0 if pd.isna(result) else float(result)

        total_bp = value("Total Interspersed Repeat", "Coverage (bp)")
        total_percent = value("Total Interspersed Repeat", "% Genome Coverage")
        genome_size = value("Total Interspersed Repeat", "Genome Size")
        if genome_size == 0:
            genome_size = float(group["Genome Size"].dropna().max())

        rows.append(
            {
                "Species": species,
                "Method": method,
                "Genome Size": genome_size,
                "Total repeat coverage (bp)": total_bp,
                "Total repeat coverage (%)": total_percent,
                "Total copy number": value(
                    "Total Interspersed Repeat", "Copy Number"
                ),
                "Unclassified coverage (bp)": value(
                    "Unclassified", "Coverage (bp)"
                ),
                "Unclassified coverage (%)": value(
                    "Unclassified", "% Genome Coverage"
                ),
                "Other-repeat coverage (bp)": value(
                    "Other (Simple Repeat, Microsatellite, RNA)",
                    "Coverage (bp)",
                ),
                "Other-repeat coverage (%)": value(
                    "Other (Simple Repeat, Microsatellite, RNA)",
                    "% Genome Coverage",
                ),
            }
        )
    return pd.DataFrame(rows)


def species_total_comparison(
    totals: pd.DataFrame,
    drayann_method: str,
) -> pd.DataFrame:
    metrics = [
        "Genome Size",
        "Total repeat coverage (bp)",
        "Total repeat coverage (%)",
        "Total copy number",
        "Unclassified coverage (bp)",
        "Unclassified coverage (%)",
        "Other-repeat coverage (bp)",
        "Other-repeat coverage (%)",
    ]
    earl = totals[totals["Method"] == "EarlGrey"][["Species", *metrics]]
    dray = totals[totals["Method"] == drayann_method][["Species", *metrics]]
    paired = earl.merge(
        dray,
        on="Species",
        how="inner",
        suffixes=(" EarlGrey", " DRayAnn"),
    )
    paired.insert(1, "DRayAnn summary", drayann_method.removeprefix("DRayAnn "))
    for metric in metrics[1:]:
        paired[f"{metric} difference (DRayAnn - EarlGrey)"] = (
            paired[f"{metric} DRayAnn"] - paired[f"{metric} EarlGrey"]
        )
        paired[f"{metric} ratio (DRayAnn / EarlGrey)"] = safe_ratio(
            paired[f"{metric} DRayAnn"], paired[f"{metric} EarlGrey"]
        )
    paired["Genome size difference (bp)"] = (
        paired["Genome Size DRayAnn"] - paired["Genome Size EarlGrey"]
    )
    return paired


def save_figure(fig: plt.Figure, outdir: Path, stem: str) -> None:
    fig.savefig(outdir / f"{stem}.png", dpi=200, bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_total_coverage(comparison: pd.DataFrame, outdir: Path, suffix: str) -> None:
    if comparison.empty:
        return
    data = comparison.sort_values("Species")
    x = list(range(len(data)))
    width = 0.38
    fig, ax = plt.subplots(figsize=(max(7, len(data) * 0.65), 5))
    ax.bar(
        [value - width / 2 for value in x],
        data["Total repeat coverage (%) EarlGrey"],
        width,
        label="EarlGrey",
        color="#4C78A8",
    )
    ax.bar(
        [value + width / 2 for value in x],
        data["Total repeat coverage (%) DRayAnn"],
        width,
        label=f"DRayAnn {suffix}",
        color="#F58518",
    )
    ax.set_xticks(x, data["Species"], rotation=45, ha="right")
    ax.set_ylabel("Total repeat coverage (% genome)")
    ax.set_title("Final annotated repeat coverage")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, outdir, f"total_repeat_coverage_{suffix}")


def plot_total_scatter(comparison: pd.DataFrame, outdir: Path, suffix: str) -> None:
    if comparison.empty:
        return
    x = comparison["Total repeat coverage (%) EarlGrey"]
    y = comparison["Total repeat coverage (%) DRayAnn"]
    lower = min(float(x.min()), float(y.min()), 0.0)
    upper = max(float(x.max()), float(y.max())) * 1.05
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, color="#6F4E7C", s=45)
    for _, row in comparison.iterrows():
        ax.annotate(
            row["Species"],
            (row["Total repeat coverage (%) EarlGrey"], row["Total repeat coverage (%) DRayAnn"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    ax.plot([lower, upper], [lower, upper], linestyle="--", color="0.5")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel("EarlGrey total repeat coverage (% genome)")
    ax.set_ylabel(f"DRayAnn {suffix} total repeat coverage (% genome)")
    ax.set_title("Total repeat coverage agreement")
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, outdir, f"total_repeat_coverage_scatter_{suffix}")


def plot_delta_heatmap(paired: pd.DataFrame, outdir: Path, suffix: str) -> None:
    class_rows = paired[~paired["TE Classification"].isin(AGGREGATE_CLASSES)]
    if class_rows.empty:
        return
    value_column = "% Genome Coverage difference (DRayAnn - EarlGrey)"
    matrix = class_rows.pivot_table(
        index="Species",
        columns="TE Classification",
        values=value_column,
        aggfunc="first",
        fill_value=0,
    )
    ordered = [column for column in CLASS_ORDER if column in matrix.columns]
    matrix = matrix.reindex(columns=ordered).sort_index()
    maximum = float(matrix.abs().to_numpy().max()) if matrix.size else 1.0
    maximum = max(maximum, 0.01)

    fig, ax = plt.subplots(
        figsize=(max(8, len(matrix.columns) * 1.15), max(4, len(matrix) * 0.45))
    )
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-maximum, vmax=maximum, aspect="auto")
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    ax.set_title(f"DRayAnn {suffix} minus EarlGrey coverage")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Difference in genome coverage (percentage points)")
    save_figure(fig, outdir, f"class_coverage_difference_heatmap_{suffix}")


def markdown_table(table: pd.DataFrame) -> str:
    def cell(value: object) -> str:
        if pd.isna(value):
            return "NA"
        if isinstance(value, float):
            value = f"{value:.6g}"
        return str(value).replace("|", "\\|").replace("\n", " ")

    header = "| " + " | ".join(cell(column) for column in table.columns) + " |"
    divider = "| " + " | ".join("---" for _ in table.columns) + " |"
    rows = [
        "| " + " | ".join(cell(value) for value in row) + " |"
        for row in table.itertuples(index=False, name=None)
    ]
    return "\n".join([header, divider, *rows])


def write_report(
    outdir: Path,
    inventory: pd.DataFrame,
    comparisons: list[pd.DataFrame],
    summary_kind: str,
) -> None:
    available = inventory[inventory["Status"] == "found"]
    species_count = available["Species"].nunique()
    lines = [
        "# EarlGrey versus DRayAnn annotation comparison",
        "",
        f"- Species discovered with at least one summary: {species_count}",
        f"- Requested DRayAnn summary: {summary_kind}",
        "- Nested rows are retained in `high_level_raw.tsv` and collapsed into their parent class in paired comparisons.",
        "- Family count is not summed when collapsing nested rows because parent and nested rows can contain the same families.",
        "",
        "This is a descriptive comparison of final annotations, not an accuracy benchmark. Without an independent truth set, higher coverage or copy number does not by itself mean better performance.",
        "",
        "## Paired total coverage",
        "",
    ]
    for comparison in comparisons:
        if comparison.empty:
            continue
        label = str(comparison["DRayAnn summary"].iloc[0])
        view = comparison[
            [
                "Species",
                "Total repeat coverage (%) EarlGrey",
                "Total repeat coverage (%) DRayAnn",
                "Total repeat coverage (%) difference (DRayAnn - EarlGrey)",
            ]
        ].copy()
        lines.extend([f"### DRayAnn {label}", "", markdown_table(view), ""])

    missing = inventory[inventory["Status"] != "found"]
    if not missing.empty:
        lines.extend(["## Missing inputs", "", markdown_table(missing), ""])

    lines.extend(
        [
            "## Output tables",
            "",
            "- `input_inventory.tsv`: discovered and missing input summaries.",
            "- `high_level_raw.tsv`: source rows, including separate nested rows.",
            "- `high_level_nested_collapsed.tsv`: parent and nested coverage combined.",
            "- `class_comparison.tsv`: class-level EarlGrey/DRayAnn differences and ratios.",
            "- `species_summary.tsv`: total, unclassified and other-repeat coverage by method.",
            "- `species_total_comparison.tsv`: paired species-level differences and ratios.",
            "",
        ]
    )
    (outdir / "comparison_report.md").write_text("\n".join(lines))


def run(args: argparse.Namespace) -> int:
    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Annotation root does not exist: {root}")
    outdir = (
        args.outdir.expanduser().resolve()
        if args.outdir is not None
        else root / "EarlGrey_vs_DRayAnn"
    )
    outdir.mkdir(parents=True, exist_ok=True)
    if not args.no_plots and plt is None:
        raise RuntimeError(
            "Plotting requires matplotlib. Install the DRayAnn project "
            "dependencies or rerun with --no-plots."
        )

    discovered = method_directories(root)
    requested_species = (
        {value.casefold() for value in args.species}
        if args.species
        else None
    )
    summary_kinds = (
        ("filtered", "unfiltered")
        if args.drayann_summary == "both"
        else (args.drayann_summary,)
    )

    inventory_rows: list[dict] = []
    tables: list[pd.DataFrame] = []
    for species_key, directories in sorted(discovered.items()):
        if requested_species is not None and species_key not in requested_species:
            continue
        species = str(directories["species"])

        earl_directory = directories.get("EarlGrey")
        earl_path = (
            find_earlgrey_summary(earl_directory, species)
            if isinstance(earl_directory, Path)
            else None
        )
        inventory_rows.append(
            {
                "Species": species,
                "Method": "EarlGrey",
                "Status": "found" if earl_path else "missing",
                "Summary file": str(earl_path or ""),
            }
        )
        if earl_path:
            tables.append(read_high_level(earl_path, species, "EarlGrey"))

        dray_directory = directories.get("DRayAnn")
        for summary_kind in summary_kinds:
            method = f"DRayAnn {summary_kind}"
            dray_path = (
                find_drayann_summary(dray_directory, species, summary_kind)
                if isinstance(dray_directory, Path)
                else None
            )
            inventory_rows.append(
                {
                    "Species": species,
                    "Method": method,
                    "Status": "found" if dray_path else "missing",
                    "Summary file": str(dray_path or ""),
                }
            )
            if dray_path:
                tables.append(read_high_level(dray_path, species, method))

    inventory = pd.DataFrame(
        inventory_rows,
        columns=("Species", "Method", "Status", "Summary file"),
    )
    inventory.to_csv(outdir / "input_inventory.tsv", sep="\t", index=False)
    if not tables:
        raise RuntimeError(f"No compatible high-level summary tables found under {root}")

    raw = pd.concat(tables, ignore_index=True, sort=False)
    raw.to_csv(outdir / "high_level_raw.tsv", sep="\t", index=False)
    collapsed = collapse_nested_rows(raw)
    collapsed.to_csv(
        outdir / "high_level_nested_collapsed.tsv", sep="\t", index=False
    )

    totals = species_summary(collapsed)
    totals.to_csv(outdir / "species_summary.tsv", sep="\t", index=False)

    class_comparisons: list[pd.DataFrame] = []
    total_comparisons: list[pd.DataFrame] = []
    for summary_kind in summary_kinds:
        method = f"DRayAnn {summary_kind}"
        class_paired = pair_with_earlgrey(collapsed, method)
        total_paired = species_total_comparison(totals, method)
        class_comparisons.append(class_paired)
        total_comparisons.append(total_paired)
        if not args.no_plots:
            plot_total_coverage(total_paired, outdir, summary_kind)
            plot_total_scatter(total_paired, outdir, summary_kind)
            plot_delta_heatmap(class_paired, outdir, summary_kind)

    pd.concat(class_comparisons, ignore_index=True, sort=False).to_csv(
        outdir / "class_comparison.tsv", sep="\t", index=False
    )
    pd.concat(total_comparisons, ignore_index=True, sort=False).to_csv(
        outdir / "species_total_comparison.tsv", sep="\t", index=False
    )
    write_report(outdir, inventory, total_comparisons, args.drayann_summary)

    paired_species = sum(len(table) for table in total_comparisons)
    print(f"Comparison written to: {outdir}")
    print(f"Paired species comparisons: {paired_species}")
    if (inventory["Status"] != "found").any():
        print(
            "Some inputs were missing; see input_inventory.tsv.",
            file=sys.stderr,
        )
    return 0


def main(argv: list[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
