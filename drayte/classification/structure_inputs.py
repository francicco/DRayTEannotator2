from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .ids import clean_family_id


def _valid_fasta(path: Path | None) -> bool:
    return bool(path and path.exists() and path.stat().st_size > 0)


def discovery_normalized_fasta(config, classification_cfg: dict | None = None) -> Path | None:
    """Return the complete normalized RepeatModeler library, if available.

    This is the broad family universe used by extension/REPAM. It should be
    preferred for structure-summary lengths over curated/final libraries, which
    may have already filtered families out.
    """
    classification_cfg = classification_cfg or {}

    explicit = classification_cfg.get("structure_base_fasta") or classification_cfg.get(
        "discovery_normalized_fasta"
    )
    if explicit:
        p = Path(explicit)
        if _valid_fasta(p):
            return p

    manifest = config.outdir_path / "discovery" / "discovery.manifest.json"
    if manifest.exists():
        try:
            data = json.loads(manifest.read_text())
            p = Path(data.get("normalized_library", ""))
            if _valid_fasta(p):
                return p
        except Exception:
            pass

    p = config.outdir_path / "discovery" / "rmodeler_dir" / f"{config.species}-families.mod.fa"
    if _valid_fasta(p):
        return p

    return None


def _records_by_clean_id(fasta: Path) -> dict[str, SeqRecord]:
    records: dict[str, SeqRecord] = {}
    for rec in SeqIO.parse(str(fasta), "fasta"):
        records[clean_family_id(rec.id)] = rec
    return records


def _iter_extension_rep_fastas(extension_final_consensuses: Path) -> Iterable[Path]:
    if not extension_final_consensuses.exists():
        return []
    return sorted(extension_final_consensuses.glob("*_rep.fa"))


def build_structure_input_fasta(
    *,
    base_fasta: Path,
    extension_final_consensuses: Path,
    output_fasta: Path,
    logger=None,
    force: bool = False,
) -> Path:
    """Build a complete FASTA for structural analyses.

    The base FASTA is usually discovery/rmodeler_dir/<species>-families.mod.fa.
    Extended consensus *_rep.fa files are preferred when present, but families
    absent from extension are retained from the base FASTA. This prevents the
    structure summary from losing consensus lengths for families that were later
    removed from curated/final annotation libraries.
    """
    if output_fasta.exists() and output_fasta.stat().st_size > 0 and not force:
        return output_fasta

    records = _records_by_clean_id(base_fasta)
    n_base = len(records)
    n_extension = 0

    for rep_fa in _iter_extension_rep_fastas(extension_final_consensuses):
        parsed = list(SeqIO.parse(str(rep_fa), "fasta"))
        if not parsed:
            continue
        rec = parsed[0]
        records[clean_family_id(rec.id)] = rec
        n_extension += 1

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    ordered = [records[k] for k in sorted(records)]
    with open(output_fasta, "w") as handle:
        SeqIO.write(ordered, handle, "fasta")

    if logger is not None:
        logger.info(
            "Wrote structure input FASTA: %s (base=%d, extension_overrides=%d, total=%d)",
            output_fasta,
            n_base,
            n_extension,
            len(records),
        )

    return output_fasta


def resolve_structure_input_fasta(
    *,
    config,
    classification_cfg: dict,
    outdir: Path,
    fallback_fasta: Path,
    logger=None,
) -> Path:
    """Resolve and, when possible, build the FASTA used for structure summary.

    This is deliberately separate from the curated library used for final
    classification. Structural evidence from extensionwork is generated for the
    broader discovery family universe.
    """
    explicit = classification_cfg.get("structure_input_fasta")
    if explicit:
        p = Path(explicit)
        if _valid_fasta(p):
            if logger is not None:
                logger.info("Using explicit structure input FASTA: %s", p)
            return p

    base = discovery_normalized_fasta(config, classification_cfg)
    if base is None:
        if logger is not None:
            logger.warning(
                "Could not find discovery normalized FASTA; using classification library for structure lengths: %s",
                fallback_fasta,
            )
        return fallback_fasta

    extension_final_consensuses = Path(
        classification_cfg.get(
            "extension_final_consensuses_dir",
            config.outdir_path / "extension" / "final_consensuses",
        )
    )

    structure_fasta = outdir / f"{config.species}.structure_input.fa"
    return build_structure_input_fasta(
        base_fasta=base,
        extension_final_consensuses=extension_final_consensuses,
        output_fasta=structure_fasta,
        logger=logger,
        force=bool(classification_cfg.get("rebuild_structure_input_fasta", False)),
    )
