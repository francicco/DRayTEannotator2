from __future__ import annotations

import json
import shutil
import tarfile
from pathlib import Path

import pandas as pd

from drayte.utils.names import safe_filename
from drayte.utils.paths import ensure_dir, stage_dir


def major_group(repeat_class: str) -> str:
    x = str(repeat_class)

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
        return "RC"
    if "Unknown" in x or "Unclassified" in x:
        return "Unclassified"

    return "Unclassified"


def link_or_copy(src: Path, dst: Path, use_symlinks: bool = True) -> bool:
    """
    Add a file to the family inspection directory.

    By default this creates symlinks to avoid duplicating large files.
    If symlinks are disabled, files are copied instead.
    """
    if not src.exists() or src.stat().st_size == 0:
        return False

    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() or dst.is_symlink():
        dst.unlink()

    if use_symlinks:
        dst.symlink_to(src.resolve())
    else:
        shutil.copy2(src, dst)

    return True


def write_metadata(path: Path, metadata: dict) -> None:
    path.write_text(json.dumps(metadata, indent=2) + "\n")


def make_family_plot(df: pd.DataFrame, family: str, outpath: Path) -> None:
    """
    Make a simple per-family divergence plot.

    This uses the refined or filtered TSV rows for one TE family.
    """
    import matplotlib.pyplot as plt

    if df.empty:
        return

    outpath.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(6, 4))
    plt.hist(df["mean_div"].dropna(), bins=30)
    plt.xlabel("RepeatMasker divergence (%)")
    plt.ylabel("Number of loci")
    plt.title(family)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def make_family_length_plot(df: pd.DataFrame, family: str, outpath: Path) -> None:
    """
    Make a simple per-family length distribution plot.
    """
    import matplotlib.pyplot as plt

    if df.empty:
        return

    outpath.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(6, 4))
    plt.hist(df["length"].dropna(), bins=30)
    plt.xlabel("Locus length (bp)")
    plt.ylabel("Number of loci")
    plt.title(family)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def archive_category(category_dir: Path, archive_path: Path) -> None:
    """
    Create a gzipped tar archive for one inspection category.

    Symlinks are dereferenced so the archive contains real files.
    """
    if not category_dir.exists():
        return

    archive_path.parent.mkdir(parents=True, exist_ok=True)

    with tarfile.open(archive_path, "w:gz", dereference=True) as tar:
        tar.add(category_dir, arcname=category_dir.name)


def find_extension_files(extension_dir: Path, safe_id: str) -> dict[str, Path]:
    """
    Locate extension-stage evidence files for one family.

    Evidence can be found in:
      - extension/final_consensuses/
      - extension/extensionwork/<family>/
      - extension/images_and_alignments/{likely_TEs,possible_SD,rejects}/

    The image/alignment folders are important because the final visual
    outputs are often stored there rather than in final_consensuses.
    """

    final_cons = extension_dir / "final_consensuses"
    work_dir = extension_dir / "extensionwork" / safe_id

    image_dirs = [
        extension_dir / "images_and_alignments" / "likely_TEs",
        extension_dir / "images_and_alignments" / "possible_SD",
        extension_dir / "images_and_alignments" / "rejects",
    ]

    candidates = {
        "rep_fa": [
            final_cons / f"{safe_id}_rep.fa",
            work_dir / f"{safe_id}_rep.fa",
        ],
        "msa_fa": [
            final_cons / f"{safe_id}_MSA_extended.fa",
            work_dir / f"{safe_id}_MSA_extended.fa",
            work_dir / "MSA-extended.fa",
            work_dir / "MSA-extended_with_rmod_cons.fa",
        ],
        "png": [
            final_cons / f"{safe_id}.png",
            work_dir / f"{safe_id}.png",
            work_dir / "img.png",
        ],
    }

    # Add direct expected paths from the final image/alignment categories.
    for d in image_dirs:
        candidates["msa_fa"].append(d / f"{safe_id}_MSA_extended.fa")
        candidates["png"].append(d / f"{safe_id}.png")

    found: dict[str, Path] = {}

    # First pass: exact expected filenames.
    for label, paths in candidates.items():
        for path in paths:
            if path.exists() and path.stat().st_size > 0:
                found[label] = path
                break

    # Second pass: glob fallback.
    # This helps with older files or partially different naming conventions.
    for d in image_dirs:
        if not d.exists():
            continue

        if "msa_fa" not in found:
            for p in d.glob(f"{safe_id}*_MSA_extended.fa"):
                if p.exists() and p.stat().st_size > 0:
                    found["msa_fa"] = p
                    break

        if "png" not in found:
            for p in d.glob(f"{safe_id}*.png"):
                if p.exists() and p.stat().st_size > 0:
                    found["png"] = p
                    break

    return found

def is_non_te_repeat(repeat_class: str, family: str = "") -> bool:
    x = f"{repeat_class} {family}".lower()
    return any(
        k in x
        for k in [
            "simple_repeat",
            "low_complexity",
            "satellite",
            "microsatellite",
            "rrna",
            "trna",
            "snrna",
            "scrna",
            "srprna",
            "srp_rna",
        ]
    )

def build_family_inspection(
    species: str,
    outdir: Path,
    refined_tsv: Path,
    curated_library: Path | None = None,
    extension_dir: Path | None = None,
    use_symlinks: bool = True,
    include_plots: bool = True,
    archive: bool = True,
    include_non_te_repeats: bool = False,
) -> dict:
    """
    Build per-family inspection directories.

    This recreates the original Ray-lab 'fordownload' philosophy, but
    with a clearer structure:

        family_inspection/
          LINE/
            familyA/
              evidence files
              plots
              metadata.json

    The refined TSV can be either:
      - annotation_refinement.tsv
      - filteredRepeats.tsv

    Use filteredRepeats.tsv for final annotation inspection.
    Use annotation_refinement.tsv for nested/overlap-aware inspection.
    """
    outdir = ensure_dir(outdir)
    manifest_rows = []

    if not refined_tsv.exists() or refined_tsv.stat().st_size == 0:
        raise FileNotFoundError(f"Missing refined TSV: {refined_tsv}")

    df = pd.read_csv(refined_tsv, sep="\t")

    required = {"family", "class", "length", "mean_div", "locus_id"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Refined TSV missing required columns: {missing}")

    extension_dir = extension_dir or Path()

    for family, fam_df in df.groupby("family", dropna=False):
        family = str(family)
        safe_id = safe_filename(family)

        repeat_class = str(fam_df["class"].iloc[0])

        if is_non_te_repeat(repeat_class, family) and not include_non_te_repeats:
            continue
	
        category = major_group(repeat_class)

        fam_dir = ensure_dir(outdir / category / safe_id)

        metadata = {
            "species": species,
            "family": family,
            "safe_family_id": safe_id,
            "class": repeat_class,
            "category": category,
            "copy_number": int(len(fam_df)),
            "total_bp": int(fam_df["length"].sum()),
            "mean_divergence": float(fam_df["mean_div"].mean()),
            "median_divergence": float(fam_df["mean_div"].median()),
        }

        # Per-family locus table
        family_tsv = fam_dir / f"{safe_id}.loci.tsv"
        fam_df.to_csv(family_tsv, sep="\t", index=False)

        # Metadata
        metadata_json = fam_dir / f"{safe_id}.metadata.json"
        write_metadata(metadata_json, metadata)

        # Extension evidence files, if available
        ext_files = find_extension_files(extension_dir, safe_id) if extension_dir.exists() else {}

        copied = []

        for label, src in ext_files.items():
            suffix = src.name.replace(f"{safe_id}", "")
            dst = fam_dir / f"{safe_id}{suffix}"

            if link_or_copy(src, dst, use_symlinks=use_symlinks):
                copied.append(label)

        # Curated library is copied globally, not split by family for now.
        # Family-specific consensus extraction can be added later.

        if include_plots:
            div_plot = fam_dir / f"{safe_id}.divergence_hist.pdf"
            len_plot = fam_dir / f"{safe_id}.length_distribution.pdf"

            make_family_plot(fam_df, family, div_plot)
            make_family_length_plot(fam_df, family, len_plot)

        manifest_rows.append(
            {
                **metadata,
                "directory": str(fam_dir),
                "evidence_files": ",".join(copied) if copied else "none",
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = outdir / f"{species}.family_inspection_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    if archive:
        archives_dir = ensure_dir(outdir / "archives")
        for category_dir in sorted([p for p in outdir.iterdir() if p.is_dir() and p.name != "archives"]):
            archive_category(
                category_dir,
                archives_dir / f"family_inspection_{category_dir.name}.tgz",
            )

    return {
        "stage": "family_inspection",
        "outdir": str(outdir),
        "manifest": str(manifest_path),
        "families": int(len(manifest)),
    }


def run(config, refinement_result: dict, logger) -> dict:
    """
    Pipeline entry point.
    """
    species = config.species
    outdir = ensure_dir(stage_dir(config.outdir_path, "family_inspection"))

    cfg = config.extra.get("family_inspection", {})

    logger.info("=" * 80)
    logger.info("STAGE: family_inspection")
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    if not cfg.get("enabled", True):
        logger.info("family_inspection disabled; skipping")
        return {
            "stage": "family_inspection",
            "enabled": False,
            "outdir": str(outdir),
        }

    # Prefer final filtered annotation if available, otherwise use overlap-aware refinement.
    refined_tsv = refinement_result.get("filtered_tsv") or refinement_result.get("refined_tsv")
    if refined_tsv is None:
        raise RuntimeError("No refined TSV available for family_inspection")

    extension_dir = config.outdir_path / "extension"

    curated_library = config.outdir_path / "curation" / "Final.RepeatModeler.Lib.fa"

    result = build_family_inspection(
        species=species,
        outdir=outdir,
        refined_tsv=Path(refined_tsv),
        curated_library=curated_library if curated_library.exists() else None,
        extension_dir=extension_dir,
        use_symlinks=bool(cfg.get("use_symlinks", True)),
        include_plots=bool(cfg.get("include_plots", True)),
        archive=bool(cfg.get("archive", True)),
        include_non_te_repeats=bool(cfg.get("include_non_te_repeats", False)),
    )

    logger.info("family_inspection completed: %s families", result["families"])

    return result
