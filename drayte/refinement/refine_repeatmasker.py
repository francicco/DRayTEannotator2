#!/usr/bin/env python3
from __future__ import annotations

import json
from statistics import mean
import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path


# ---------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------

@dataclass
class RMHit:
    """
    One raw RepeatMasker hit parsed from a RepeatMasker .out file.

    This represents one alignment fragment, not necessarily one biological
    TE insertion. RepeatMasker often reports a single TE insertion as several
    adjacent fragments. The goal of the refinement step is to group compatible
    fragments into larger TE loci.
    """

    # RepeatMasker Smith-Waterman score
    sw_score: int

    # Percent divergence reported by RepeatMasker
    perc_div: float

    # Genome/scaffold/contig name
    query: str

    # Genomic coordinates, normalized so q_start <= q_end
    q_start: int
    q_end: int

    # Strand: '+' or '-'
    # RepeatMasker reports 'C' for complement; we convert that to '-'
    strand: str

    # Repeat family name, e.g. rnd-1_family-10
    repeat_name: str

    # Repeat class/family, e.g. LINE/RTE, LTR/Gypsy, Unknown
    repeat_class: str

    # Coordinates on the repeat consensus.
    # Kept as strings because RepeatMasker can use parentheses, e.g. "(123)"
    r_start: str
    r_end: str
    r_left: str

    # RepeatMasker hit ID, if present
    rm_id: str

    # Length of the genomic hit
    length: int

    # Whether this hit is almost entirely contained inside a larger hit
    nested: bool = False

    # Stable internal index based on input order
    raw_index: int = 0


@dataclass
class TELocus:
    """
    One refined TE locus.

    A TELocus can consist of one RepeatMasker hit or several merged
    compatible hits. This is closer to a biological insertion than the raw
    RepeatMasker fragments.
    """

    locus_id: str
    query: str
    start: int
    end: int
    strand: str
    repeat_name: str
    repeat_class: str
    hits: list[RMHit] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start + 1


# ---------------------------------------------------------------------
# Parsing RepeatMasker output
# ---------------------------------------------------------------------

def parse_repeatmasker_out(path: Path) -> list[RMHit]:
    """
    Parse a standard RepeatMasker .out file.

    RepeatMasker .out has a header followed by rows like:

    SW perc perc perc query begin end (left) strand repeat class/family ...

    This parser skips header lines and malformed rows. It converts each
    valid line into an RMHit object.

    Important:
    - It parses the .out file directly rather than the GFF because .out
      contains more repeat-coordinate information.
    - It normalizes genomic coordinates.
    - It converts RepeatMasker's 'C' strand to '-'.
    """

    hits: list[RMHit] = []

    with path.open() as handle:
        for line in handle:
            # Skip empty lines
            if not line.strip():
                continue

            # Skip RepeatMasker header lines
            if line.lstrip().startswith(("SW", "score", "---")):
                continue

            parts = line.strip().split()

            # Standard RepeatMasker rows should have at least 14 columns.
            # Some have a 15th ID column.
            if len(parts) < 14:
                continue

            try:
                sw_score = int(parts[0])
                perc_div = float(parts[1])
                query = parts[4]
                q_start = int(parts[5])
                q_end = int(parts[6])
                strand = parts[8]
                repeat_name = parts[9]
                repeat_class = parts[10]
                r_start = parts[11]
                r_end = parts[12]
                r_left = parts[13]

                # RepeatMasker ID column, if present.
                # If absent, create a synthetic ID to avoind redundancies.
                rm_id = f"RMhit_{len(hits)+1}"

            except ValueError:
                # Skip rows that cannot be parsed numerically.
                continue

            # RepeatMasker uses C for complement strand.
            if strand == "C":
                strand = "-"

            # Store hit with normalized genome coordinates.
            hits.append(
                RMHit(
                    sw_score=sw_score,
                    perc_div=perc_div,
                    query=query,
                    q_start=min(q_start, q_end),
                    q_end=max(q_start, q_end),
                    strand=strand,
                    repeat_name=repeat_name,
                    repeat_class=repeat_class,
                    r_start=r_start,
                    r_end=r_end,
                    r_left=r_left,
                    rm_id=rm_id,
                    length=abs(q_end - q_start) + 1,
                    raw_index=len(hits) + 1,
                )
            )

    if not hits:
        raise ValueError(f"No RepeatMasker hits parsed from {path}")

    return hits


# ---------------------------------------------------------------------
# Nested hit detection
# ---------------------------------------------------------------------

def mark_nested_hits(hits: list[RMHit], min_overlap_fraction: float = 0.95) -> None:
    """
    Mark hits that are almost completely contained within larger hits.

    A hit is marked as nested if:
    - another hit on the same scaffold fully contains it
    - the containing hit is longer
    - the contained hit is covered by at least min_overlap_fraction

    This is deliberately simple. It catches cases where RepeatMasker reports
    a small hit entirely inside a larger TE annotation.

    By default, nested hits are excluded from the refined main GFF but written
    separately to nestedRepeats.gff3.
    """

    by_query: dict[str, list[RMHit]] = {}

    # Group hits by scaffold/contig
    for h in hits:
        by_query.setdefault(h.query, []).append(h)

    for qhits in by_query.values():
        # Sort so larger intervals starting at same coordinate come first
        qhits.sort(key=lambda h: (h.q_start, -(h.q_end - h.q_start)))

        for i, h in enumerate(qhits):
            for j, other in enumerate(qhits):
                if i == j:
                    continue

                # Check whether "other" fully contains "h"
                if other.q_start <= h.q_start and other.q_end >= h.q_end:
                    if other.length <= h.length:
                        continue

                    # Since h is fully contained, overlap == h.length.
                    overlap = h.length

                    if overlap / h.length >= min_overlap_fraction:
                        h.nested = True
                        break


# ---------------------------------------------------------------------
# Fragment merging
# ---------------------------------------------------------------------

def compatible_for_merge(a: RMHit, b: RMHit, max_gap: int, max_consensus_overlap=50, max_consensus_jump=10000) -> bool:
    """
    Decide whether two adjacent RepeatMasker hits should be merged.

    Conservative merge criteria:
    - same scaffold
    - same family name
    - same class
    - same strand
    - non-overlapping
    - distance between them <= max_gap

    This intentionally avoids merging:
    - different families
    - different strands
    - overlapping hits
    - distant fragments

    The current version does not yet use repeat-coordinate compatibility
    (r_start/r_end/r_left). That could be added later to avoid merging
    biologically inconsistent fragments.
    """

    if a.query != b.query:
        return False

    if a.repeat_name != b.repeat_name:
        return False

    if a.repeat_class != b.repeat_class:
        return False

    if a.strand != b.strand:
        return False

    # Distance between two genomic fragments.
    # Example:
    # a: 100-200
    # b: 251-300
    # gap = 50
    gap = b.q_start - a.q_end - 1

    # Negative gap means overlap; current conservative mode does not merge overlaps.
    if gap < 0:
        return False

    if gap > max_gap:
        return False

    return True

def is_ignorable_intervening_hit(h: RMHit) -> bool:
    """
    Decide whether a hit can be ignored when merging TE fragments in 'loose' mode.

    Concept:
    In loose mode, we allow merging two fragments of the same TE even if
    something lies between them — BUT only if that something is biologically
    uninformative or likely spurious.

    Typical "ignorable" elements:
    - simple repeats (e.g. (AT)n)
    - low complexity regions
    - satellites / microsatellites
    - small RNA annotations (tRNA, rRNA, etc.)
    - already-detected nested fragments

    These elements often interrupt TE annotations but are not true boundaries.

    Implementation details:
    - We inspect both RepeatMasker class (repeat_class)
      and repeat name (repeat_name), because annotation conventions vary.
    - Matching is done in lowercase to avoid case issues.

    Limitations:
    - This is heuristic and string-based (like RepeatCraft).
    - Could be refined later with a controlled vocabulary.
    """

    cls = h.repeat_class.lower()
    name = h.repeat_name.lower()

    return (
        # Already marked nested → safe to ignore
        h.nested

        # Simple repeats / low complexity
        or "simple_repeat" in cls
        or "simple" in cls
        or "low_complexity" in cls

        # Satellites / microsatellites
        or "satellite" in cls
        or "microsatellite" in cls

        # RNA-derived annotations (non-TE interruptions)
        or "trna" in cls
        or "rrna" in cls
        or "snrna" in cls
        or "scrna" in cls
        or "srp_rna" in cls

        # Fallback checks on repeat name
        or "simple_repeat" in name
        or "low_complexity" in name
    )

def compatible_for_loose_merge(
    current: TELocus,
    candidate: RMHit,
    intervening: list[RMHit],
    max_gap: int,
    max_consensus_overlap=50,
    max_consensus_jump=10000
) -> bool:
    """
    Decide whether a new hit can be merged into the current TE locus
    under 'loose' mode.

    Difference from strict mode:
    - strict: only merge immediately adjacent compatible fragments
    - loose: allow intervening hits IF they are ignorable

    Inputs:
    - current: the TE locus being built
    - candidate: the next hit we want to evaluate
    - intervening: list of hits encountered since the last accepted fragment
    - max_gap: maximum allowed genomic distance

    Logic:

    Step 1: Ensure same biological identity
        - same scaffold (query)
        - same repeat family (repeat_name)
        - same class
        - same strand

    Step 2: Check genomic proximity
        - distance between current locus end and candidate start must be:
            0 ≤ gap ≤ max_gap

    Step 3: Validate intervening elements
        - ALL intervening hits must be ignorable
        - if even one is not ignorable → do NOT merge

    Important:
    - This prevents merging across unrelated TEs.
    - This is the key difference vs strict mode.

    Limitations:
    - Does not yet check consistency of repeat coordinates (r_start/r_end).
    """

    # Compare against last fragment in current locus
    last = current.hits[-1]

    # --- Identity checks ---
    if last.query != candidate.query:
        return False

    if last.repeat_name != candidate.repeat_name:
        return False

    if last.repeat_class != candidate.repeat_class:
        return False

    if last.strand != candidate.strand:
        return False

    # --- Distance check ---
    gap = candidate.q_start - current.end - 1

    if gap < 0:
        # Overlapping fragments are not merged in this version
        return False

    if gap > max_gap:
        return False

    # --- Intervening check ---
    # Only allow merge if ALL intervening hits are ignorable
    return all(is_ignorable_intervening_hit(h) for h in intervening)

def build_loci(
    hits: list[RMHit],
    species: str,
    max_gap: int = 150,
    include_nested: bool = False,
    mode: str = "strict",
    max_consensus_overlap: int = 50,
    max_consensus_jump: int = 10000,
) -> list[TELocus]:
    """
    Build TE loci by merging RepeatMasker hits.

    Two modes:
    - strict: merge only directly adjacent compatible fragments
    - loose: merge across ignorable intervening hits

    Workflow:

    1. Sort hits by genomic coordinates
    2. Iterate linearly
    3. Maintain:
        - current locus
        - buffer of intervening hits
    4. Decide whether to:
        - extend current locus
        - buffer hit
        - start new locus
    """

    if mode not in {"strict", "loose"}:
        raise ValueError("mode must be either 'strict' or 'loose'")

    # Sort by genome coordinate
    sorted_hits = sorted(
        hits,
        key=lambda h: (h.query, h.q_start, h.q_end, h.repeat_name),
    )

    loci: list[TELocus] = []
    current: TELocus | None = None

    # Holds hits encountered between mergeable fragments
    intervening_buffer: list[RMHit] = []

    for h in sorted_hits:

        # ----------------------------------------------------------
        # Handle nested hits
        # ----------------------------------------------------------
        if h.nested and not include_nested:
            # Do not include in main loci
            # but keep track of them as intervening candidates
            intervening_buffer.append(h)
            continue

        # ----------------------------------------------------------
        # Initialize first locus
        # ----------------------------------------------------------
        if current is None:

            # In loose mode, skip leading ignorable hits
            if mode == "loose" and is_ignorable_intervening_hit(h) and not include_nested:
                intervening_buffer.append(h)
                continue

            current = TELocus(
                locus_id="",
                query=h.query,
                start=h.q_start,
                end=h.q_end,
                strand=h.strand,
                repeat_name=h.repeat_name,
                repeat_class=h.repeat_class,
                hits=[h],
            )

            intervening_buffer = []
            continue

        # ----------------------------------------------------------
        # STRICT MODE
        # ----------------------------------------------------------
        if mode == "strict":

            if compatible_for_merge(current.hits[-1], h, max_gap, max_consensus_overlap=max_consensus_overlap, max_consensus_jump=max_consensus_jump):
                # Direct extension
                current.end = max(current.end, h.q_end)
                current.hits.append(h)
            else:
                # Break locus
                loci.append(current)

                current = TELocus(
                    locus_id="",
                    query=h.query,
                    start=h.q_start,
                    end=h.q_end,
                    strand=h.strand,
                    repeat_name=h.repeat_name,
                    repeat_class=h.repeat_class,
                    hits=[h],
                )

        # ----------------------------------------------------------
        # LOOSE MODE
        # ----------------------------------------------------------
        else:

            # Case 1: direct compatibility (same as strict)
            if compatible_for_merge(current.hits[-1], h, max_gap, max_consensus_overlap=max_consensus_overlap, max_consensus_jump=max_consensus_jump):
                current.end = max(current.end, h.q_end)
                current.hits.append(h)

                # Reset buffer after successful merge
                intervening_buffer = []

            # Case 2: merge across intervening ignorable hits
            elif compatible_for_loose_merge(current, h, intervening_buffer, max_gap, max_consensus_overlap=max_consensus_overlap, max_consensus_jump=max_consensus_jump):
                current.end = max(current.end, h.q_end)
                current.hits.append(h)

                # Clear buffer once used
                intervening_buffer = []

            # Case 3: hit is ignorable → buffer it
            elif is_ignorable_intervening_hit(h) and h.query == current.query:
                intervening_buffer.append(h)

            # Case 4: real boundary → start new locus
            else:
                loci.append(current)

                current = TELocus(
                    locus_id="",
                    query=h.query,
                    start=h.q_start,
                    end=h.q_end,
                    strand=h.strand,
                    repeat_name=h.repeat_name,
                    repeat_class=h.repeat_class,
                    hits=[h],
                )

                intervening_buffer = []

    # Final locus
    if current is not None:
        loci.append(current)

    # ----------------------------------------------------------
    # Assign stable IDs
    # ----------------------------------------------------------
    width = max(7, len(str(len(loci))))
    for i, locus in enumerate(loci, start=1):
        locus.locus_id = f"{species}_TE{str(i).zfill(width)}"

    return loci

def mark_locus_overlaps(loci: list[TELocus]) -> dict[str, list[str]]:
    overlaps: dict[str, list[str]] = {l.locus_id: [] for l in loci}

    by_seqid: dict[str, list[TELocus]] = {}
    for locus in loci:
        by_seqid.setdefault(locus.query, []).append(locus)

    for seq_loci in by_seqid.values():
        seq_loci.sort(key=lambda x: (x.start, x.end))

        for i, a in enumerate(seq_loci):
            for b in seq_loci[i + 1:]:
                if b.start > a.end:
                    break
                if a.start <= b.end and b.start <= a.end:
                    overlaps[a.locus_id].append(b.locus_id)
                    overlaps[b.locus_id].append(a.locus_id)

    return overlaps

# ---------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------

def gff_escape(value: str) -> str:
    """
    Escape problematic characters in GFF3 attributes.

    GFF3 attributes use semicolon-separated key=value pairs.
    Characters such as ';', '=', '&', ',' and spaces can break parsers.
    """

    return (
        str(value)
        .replace(";", "%3B")
        .replace("=", "%3D")
        .replace("&", "%26")
        .replace(",", "%2C")
        .replace(" ", "_")
    )


def write_gff3(loci: list[TELocus], path: Path, overlaps: dict[str, list[str]] | None = None) -> None:
    """
    Write refined TE loci as GFF3.

    Each locus is written as one repeat_region feature.
    The attributes include:
    - ID: stable locus ID
    - Name: repeat family name
    - Family: repeat family name
    - Class: RepeatMasker class/family
    - Fragments: number of raw RepeatMasker hits merged
    - MeanDiv: mean RepeatMasker divergence of fragments
    - RawIDs: original RepeatMasker hit IDs
    """

    with path.open("w") as out:
        out.write("##gff-version 3\n")

        for locus in loci:
            overlap_ids = overlaps.get(locus.locus_id, []) if overlaps else []

            attrs = {
                "ID": locus.locus_id,
                "Name": locus.repeat_name,
                "Family": locus.repeat_name,
                "Class": locus.repeat_class,
                "Fragments": str(len(locus.hits)),
                "MeanDiv": f"{sum(h.perc_div for h in locus.hits) / len(locus.hits):.3f}",
                "RawIDs": ",".join(h.rm_id for h in locus.hits),
                "Overlaps": "Yes" if overlap_ids else "No",
                "OverlapIDs": ",".join(overlap_ids) if overlap_ids else "None",
            }

            attr_str = ";".join(f"{k}={gff_escape(v)}" for k, v in attrs.items())

            out.write(
                "\t".join(
                    [
                        locus.query,
                        "DRayTE_refine",
                        "repeat_region",
                        str(locus.start),
                        str(locus.end),
                        ".",
                        locus.strand,
                        ".",
                        attr_str,
                    ]
                )
                + "\n"
            )

def write_refinement_stats(
    loci: list[TELocus],
    hits: list[RMHit],
    path: Path,
    params: dict,
    overlaps: dict[str, list[str]],
) -> dict:
    fragment_counts = [len(l.hits) for l in loci]
    overlapping_loci = sum(1 for locus_id, ids in overlaps.items() if ids)

    stats = {
        **params,
        "raw_hits": len(hits),
        "nested_hits": sum(1 for h in hits if h.nested),
        "refined_loci": len(loci),
        "single_fragment_loci": sum(1 for x in fragment_counts if x == 1),
        "merged_loci": sum(1 for x in fragment_counts if x > 1),
        "max_fragments_per_locus": max(fragment_counts) if fragment_counts else 0,
        "mean_fragments_per_locus": round(mean(fragment_counts), 4) if fragment_counts else 0,
        "total_refined_bp": sum(l.length for l in loci),
        "overlapping_loci": overlapping_loci,
    }

    with path.open("w") as out:
        out.write("metric\tvalue\n")
        for k, v in stats.items():
            out.write(f"{k}\t{v}\n")

    return stats

def write_manifest(path: Path, result: dict) -> None:
    path.write_text(json.dumps(result, indent=2) + "\n")

def write_bed(loci: list[TELocus], path: Path) -> None:
    """
    Write refined TE loci as BED6.

    BED is 0-based half-open, so GFF start is converted:
    BED start = GFF start - 1
    BED end   = GFF end
    """

    with path.open("w") as out:
        for locus in loci:
            out.write(
                "\t".join(
                    [
                        locus.query,
                        str(locus.start - 1),
                        str(locus.end),
                        locus.locus_id,
                        "0",
                        locus.strand,
                    ]
                )
                + "\n"
            )


def major_class_from_repeat_class(repeat_class: str) -> str:
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
    if any(k in x for k in ["Simple_repeat", "Low_complexity", "Satellite", "RNA", "rRNA", "tRNA", "snRNA", "srpRNA"]):
        return "Other (Simple Repeat, Microsatellite, RNA)"
    if "Unknown" in x or "Unclassified" in x:
        return "Unclassified"
    return "Unclassified"

def write_tsv(
    loci: list[TELocus],
    path: Path,
    overlaps: dict[str, list[str]] | None = None,
) -> None:
    """
    Write a tabular summary of refined loci.

    This is useful for debugging and downstream statistics.
    """

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow(
            [
                "locus_id",
                "seqid",
                "start",
                "end",
                "length",
                "strand",
                "family",
                "class",
                "classification",
                "n_fragments",
                "mean_div",
                "raw_ids",
                "overlaps",
                "overlap_ids",
            ]
        )

        for locus in loci:
            overlap_ids = overlaps.get(locus.locus_id, []) if overlaps else []
            has_overlap = bool(overlap_ids)

            base_class = major_class_from_repeat_class(locus.repeat_class)
            classification = f"{base_class}-nested" if has_overlap else base_class

            writer.writerow(
                [
                    locus.locus_id,
                    locus.query,
                    locus.start,
                    locus.end,
                    locus.length,
                    locus.strand,
                    locus.repeat_name,
                    locus.repeat_class,
                    classification,
                    len(locus.hits),
                    round(sum(h.perc_div for h in locus.hits) / len(locus.hits), 3),
                    ",".join(h.rm_id for h in locus.hits),
                    "Yes" if has_overlap else "No",
                    ",".join(overlap_ids) if overlap_ids else "None",
                ]
            )

def write_nested_gff(hits: list[RMHit], path: Path, species: str) -> None:
    """
    Write hits marked as nested to a separate GFF3 file.

    These are excluded from the main refined annotation by default, but
    preserved for inspection.
    """

    nested = [h for h in hits if h.nested]

    with path.open("w") as out:
        out.write("##gff-version 3\n")

        for i, h in enumerate(nested, start=1):
            attrs = {
                "ID": f"{species}_nested_{i}",
                "Name": h.repeat_name,
                "Family": h.repeat_name,
                "Class": h.repeat_class,
                "RawID": h.rm_id,
            }

            attr_str = ";".join(f"{k}={gff_escape(v)}" for k, v in attrs.items())

            out.write(
                "\t".join(
                    [
                        h.query,
                        "DRayTE_refine",
                        "nested_repeat",
                        str(h.q_start),
                        str(h.q_end),
                        str(h.sw_score),
                        h.strand,
                        ".",
                        attr_str,
                    ]
                )
                + "\n"
            )


# ---------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------

def validate_gff3(path: Path) -> None:
    """
    Fail-fast validation for refined GFF3.

    This prevents the EarlGrey-style problem where downstream steps continue
    even if attributes are missing or malformed.
    """

    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty GFF3: {path}")

    total = 0
    bad = 0

    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue

            total += 1
            parts = line.rstrip("\n").split("\t")

            if len(parts) != 9:
                bad += 1
                continue

            if "ID=" not in parts[8] or "Name=" not in parts[8]:
                bad += 1

    if total == 0:
        raise ValueError(f"No feature rows in GFF3: {path}")

    if bad:
        raise ValueError(f"Malformed GFF3 rows in {path}: {bad}/{total}")


# ---------------------------------------------------------------------
# Main refinement function
# ---------------------------------------------------------------------

def refine_repeatmasker(
    rmout: Path,
    species: str,
    outdir: Path,
    max_gap: int = 150,
    include_nested: bool = False,
    mode: str = "strict",
    max_consensus_overlap: int = 50,
    max_consensus_jump: int = 10000,
) -> dict:
    """
    Run the full refinement workflow.

    Input:
      - RepeatMasker .out
      - species prefix
      - output directory

    Output:
      - refined GFF3
      - refined BED
      - refinement TSV
      - nested-repeat GFF3

    This is the function used both by:
      - the standalone TE-refine command
      - the DRayTE pipeline stage
    """

    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Parse raw RepeatMasker hits
    hits = parse_repeatmasker_out(rmout)

    # 2. Mark fully nested hits
    mark_nested_hits(hits)

    # 3. Build refined loci by conservative fragment merging
    loci = build_loci(
        hits=hits,
        species=species,
        max_gap=max_gap,
        include_nested=include_nested,
        mode=mode,
        max_consensus_overlap=max_consensus_overlap,
        max_consensus_jump=max_consensus_jump
    )

    overlaps = mark_locus_overlaps(loci)

    stats_tsv = outdir / f"{species}.annotation_refinement.stats.tsv"
    manifest_json = outdir / f"{species}.annotation_refinement.manifest.json"

    # 4. Define outputs
    refined_gff = outdir / f"{species}.refinedRepeats.gff3"
    refined_bed = outdir / f"{species}.refinedRepeats.bed"
    refined_tsv = outdir / f"{species}.annotation_refinement.tsv"
    nested_gff = outdir / f"{species}.nestedRepeats.gff3"

    # 5. Write outputs
    write_gff3(loci, refined_gff, overlaps=overlaps)
    write_bed(loci, refined_bed)
    write_tsv(loci, refined_tsv, overlaps=overlaps)
    write_nested_gff(hits, nested_gff, species)

    # 6. Validate final GFF3
    validate_gff3(refined_gff)

	# 7. Write write refinement stats
    params = {
        "mode": mode,
        "max_gap": max_gap,
        "include_nested": include_nested,
        "max_consensus_overlap": max_consensus_overlap,
        "max_consensus_jump": max_consensus_jump,
    }
    
    stats = write_refinement_stats(
        loci=loci,
        hits=hits,
        path=stats_tsv,
        params=params,
        overlaps=overlaps,
    )

    # 8. Return metadata for manifest/reporting
    result = {
        **params,
        **stats,
        "refined_gff": str(refined_gff),
        "refined_bed": str(refined_bed),
        "refined_tsv": str(refined_tsv),
        "nested_gff": str(nested_gff),
        "stats_tsv": str(stats_tsv),
        "manifest_json": str(manifest_json),
    }
    
    write_manifest(manifest_json, result)
    
    return result


# ---------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------

def main() -> None:
    """
    Standalone CLI.

    Example:
      TE-refine \
        --rmout Cant.genome.v1.fasta.out \
        --species Cant \
        --outdir Cant.annotation_refinement \
        --max-gap 150
    """

    parser = argparse.ArgumentParser(
        description=(
            "DRayTE RepeatMasker annotation refinement: "
            "conservative or loose RepeatCraft-like defragmentation."
        )
    )

    parser.add_argument(
        "--rmout",
        required=True,
        type=Path,
        help="RepeatMasker .out file",
    )

    parser.add_argument(
        "--species",
        required=True,
        help="Species/sample prefix",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory",
    )

    parser.add_argument(
        "--max-gap",
        type=int,
        default=150,
        help="Maximum gap between same-family fragments to merge. Default: 150",
    )

    parser.add_argument(
        "--include-nested",
        action="store_true",
        help="Keep nested hits in refined GFF instead of excluding them",
    )

    parser.add_argument(
        "--mode",
        choices=["strict", "loose"],
        default="strict",
        help="Fragment merging mode. strict = adjacent only; loose = merge across ignorable intervening hits. Default: strict",
    )

    parser.add_argument(
        "--max-consensus-overlap",
        type=int,
        default=50,
        help="Maximum allowed overlap between repeat-consensus coordinates when merging fragments. Default: 50",
    )
    
    parser.add_argument(
        "--max-consensus-jump",
        type=int,
        default=10000,
        help="Maximum allowed jump between repeat-consensus coordinates when merging fragments. Default: 10000",
    )

    args = parser.parse_args()

    result = refine_repeatmasker(
        rmout=args.rmout,
        species=args.species,
        outdir=args.outdir,
        max_gap=args.max_gap,
        include_nested=args.include_nested,
		mode=args.mode,
        max_consensus_overlap=args.max_consensus_overlap,
        max_consensus_jump=args.max_consensus_jump,
    )

    # Print summary in machine-readable key/value format
    for k, v in result.items():
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
