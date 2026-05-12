from __future__ import annotations

import csv
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from .ids import clean_family_id


@dataclass
class RepamRange:
    seq_id: str
    chrom: str
    start: int
    end: int
    strand: str
    extended_left: int
    extended_right: int
    anchor_range: str = ""


@dataclass
class CopyTSDHit:
    seq_id: str
    tsd_len: int
    tsd_seq: str
    left_pos: int
    right_pos: int
    left_offset: int = 0
    right_offset: int = 0
    boundary_distance: int = 0


@dataclass
class RepamTSDCall:
    family_id: str
    tsd_present: bool
    tsd_len: int
    tsd_seq: str
    tsd_support: float
    n_copies_checked: int
    n_copies_with_tsd: int
    n_distinct_tsd_hits: int
    source: str = "repam"

    # TSD support is split into length-level and exact-sequence support.
    # The public tsd_support field is kept for backwards compatibility and
    # stores the biologically more useful length-level support.
    tsd_len_support: float = 0.0
    tsd_exact_support: float = 0.0
    tsd_confidence: str = "NONE"
    n_copies_supporting_len: int = 0
    n_copies_supporting_exact: int = 0


def is_low_complexity(
    seq: str,
    max_single_base_fraction: float = 0.80,
) -> bool:
    seq = seq.upper()
    seq = "".join(x for x in seq if x in {"A", "C", "G", "T"})

    if not seq:
        return True

    if len(set(seq)) < 2:
        return True

    max_base_fraction = max(
        seq.count(base) / len(seq)
        for base in "ACGT"
    )

    return max_base_fraction >= max_single_base_fraction


def parse_repam_metadata(meta: str) -> dict[str, str]:
    values: dict[str, str] = {}

    for item in meta.split(","):
        item = item.strip()

        if "=" not in item:
            continue

        key, value = item.split("=", 1)
        values[key.strip()] = value.strip()

    return values


def parse_repam_ranges(path: str | Path) -> list[RepamRange]:
    """
    Parse repam-ranges.tsv.

    Expected rows resemble:

    chrom    start    end    strand    n=0,anchor_range=...,extended_left=...,extended_right=...

    The n= index is used to associate the range with the corresponding
    repam-repseq FASTA record if needed.
    """

    path = Path(path)
    ranges: list[RepamRange] = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) < 5:
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            strand = parts[3]
            meta = parts[4]

            md = parse_repam_metadata(meta)

            n = md.get("n", str(len(ranges)))
            seq_id = f"n={n}"

            extended_left = int(md.get("extended_left", 0))
            extended_right = int(md.get("extended_right", 0))

            ranges.append(
                RepamRange(
                    seq_id=seq_id,
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    extended_left=extended_left,
                    extended_right=extended_right,
                    anchor_range=md.get("anchor_range", ""),
                )
            )

    return ranges


def fasta_records_by_index(path: str | Path) -> dict[str, str]:
    """
    Return FASTA records keyed by n=<index> and by raw record id.

    REPAM output usually keeps records in the same order as repam-ranges.tsv.
    This function therefore provides both:
      n=0, n=1, ...
      actual FASTA record ID
    """

    records: dict[str, str] = {}

    for i, rec in enumerate(SeqIO.parse(str(path), "fasta")):
        seq = str(rec.seq).upper()

        records[f"n={i}"] = seq
        records[rec.id] = seq

    return records


def best_tsd_for_copy(
    seq: str,
    extended_left: int,
    extended_right: int,
    min_len: int = 2,
    max_len: int = 15,
    boundary_slop: int = 3,
) -> CopyTSDHit | None:
    """
    Infer a TSD for one REPAM-extended copy.

    Sequence layout is assumed to be approximately:

        left_extension | TE_body | right_extension

    In practice, consensus boundaries are often shifted by a few bases.
    Therefore the search allows a small independent offset around both
    predicted boundaries. The best hit is the longest exact non-low-complexity
    repeat with the smallest total boundary displacement.
    """

    seq = seq.upper()

    if extended_left <= 0 or extended_right <= 0:
        return None

    te_start = extended_left
    te_end = len(seq) - extended_right

    if te_start <= 0 or te_end <= te_start:
        return None

    candidates: list[CopyTSDHit] = []

    for k in range(max_len, min_len - 1, -1):
        for left_offset in range(-boundary_slop, boundary_slop + 1):
            left_boundary = te_start + left_offset

            if left_boundary < k:
                continue

            for right_offset in range(-boundary_slop, boundary_slop + 1):
                right_boundary = te_end + right_offset

                if right_boundary < 0:
                    continue

                if right_boundary + k > len(seq):
                    continue

                left = seq[left_boundary - k:left_boundary]
                right = seq[right_boundary:right_boundary + k]

                if "N" in left or "N" in right:
                    continue

                if left != right:
                    continue

                if is_low_complexity(left):
                    continue

                candidates.append(
                    CopyTSDHit(
                        seq_id="",
                        tsd_len=k,
                        tsd_seq=left,
                        left_pos=left_boundary - k + 1,
                        right_pos=right_boundary + 1,
                        left_offset=left_offset,
                        right_offset=right_offset,
                        boundary_distance=abs(left_offset) + abs(right_offset),
                    )
                )


    if not candidates:
        return None

    return sorted(
        candidates,
        key=lambda x: (
            x.boundary_distance,
            abs(x.left_offset - x.right_offset),
            -x.tsd_len,
            x.tsd_seq,
        ),
    )[0]


def infer_family_tsd_from_repam(
    repam_ranges: str | Path,
    repam_repseq: str | Path,
    family_id: str | None = None,
    min_len: int = 2,
    max_len: int = 15,
    min_copies: int = 3,
    min_support: float = 0.30,
    boundary_slop: int = 3,
) -> RepamTSDCall:
    repam_ranges = Path(repam_ranges)
    repam_repseq = Path(repam_repseq)

    if family_id is None:
        family_id = clean_family_id(repam_ranges.parent.name)

    ranges = parse_repam_ranges(repam_ranges)
    seqs = fasta_records_by_index(repam_repseq)

    observed: list[tuple[int, str]] = []
    n_copies_checked = 0

    for rg in ranges:
        seq = seqs.get(rg.seq_id)

        if seq is None:
            continue

        n_copies_checked += 1

        hit = best_tsd_for_copy(
            seq=seq,
            extended_left=rg.extended_left,
            extended_right=rg.extended_right,
            min_len=min_len,
            max_len=max_len,
            boundary_slop=boundary_slop,
        )

        if hit is None:
            continue

        hit.seq_id = rg.seq_id
        observed.append((hit.tsd_len, hit.tsd_seq))

    if n_copies_checked < min_copies or not observed:
        return RepamTSDCall(
            family_id=family_id,
            tsd_present=False,
            tsd_len=0,
            tsd_seq="",
            tsd_support=0.0,
            n_copies_checked=n_copies_checked,
            n_copies_with_tsd=len(observed),
            n_distinct_tsd_hits=len(set(observed)),
            tsd_len_support=0.0,
            tsd_exact_support=0.0,
            tsd_confidence="NONE",
            n_copies_supporting_len=0,
            n_copies_supporting_exact=0,
        )

    exact_counts = Counter(observed)
    len_counts = Counter(length for length, _seq in observed)

    best_len, best_len_n = len_counts.most_common(1)[0]
    seq_counts_for_best_len = Counter(
        seq for length, seq in observed
        if length == best_len
    )
    best_seq, best_exact_n = seq_counts_for_best_len.most_common(1)[0]

    len_support = best_len_n / n_copies_checked
    exact_support = best_exact_n / n_copies_checked

    tsd_present = len_support >= min_support

    if not tsd_present:
        confidence = "NONE"
    elif exact_support >= 0.50 or (len_support >= 0.70 and exact_support >= 0.30):
        confidence = "HIGH"
    elif len_support >= 0.50 or exact_support >= 0.30:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"

    return RepamTSDCall(
        family_id=family_id,
        tsd_present=tsd_present,
        tsd_len=best_len,
        tsd_seq=best_seq,
        # Keep tsd_support as the length-level support. TSD length is often
        # more useful than exact sequence for TIR superfamily inference.
        tsd_support=round(len_support, 3),
        n_copies_checked=n_copies_checked,
        n_copies_with_tsd=len(observed),
        n_distinct_tsd_hits=len(exact_counts),
        tsd_len_support=round(len_support, 3),
        tsd_exact_support=round(exact_support, 3),
        tsd_confidence=confidence,
        n_copies_supporting_len=best_len_n,
        n_copies_supporting_exact=best_exact_n,
    )


def infer_tsds_from_extensionwork(
    extensionwork_dir: str | Path,
    min_len: int = 2,
    max_len: int = 15,
    min_copies: int = 3,
    min_support: float = 0.30,
    boundary_slop: int = 3,
) -> list[RepamTSDCall]:
    extensionwork_dir = Path(extensionwork_dir)

    calls: list[RepamTSDCall] = []

    for family_dir in sorted(extensionwork_dir.iterdir()):
        if not family_dir.is_dir():
            continue

        repam_ranges = family_dir / "repam-ranges.tsv"
        repam_repseq = family_dir / "repam-repseq.fa"

        if not repam_ranges.exists():
            continue

        if not repam_repseq.exists():
            continue

        calls.append(
            infer_family_tsd_from_repam(
                repam_ranges=repam_ranges,
                repam_repseq=repam_repseq,
                family_id=clean_family_id(family_dir.name),
                min_len=min_len,
                max_len=max_len,
                min_copies=min_copies,
                min_support=min_support,
                boundary_slop=boundary_slop,
            )
        )

    return calls


def write_repam_tsd_tsv(
    calls: list[RepamTSDCall],
    outfile: str | Path,
) -> None:
    fields = [
        "family_id",
        "tsd_present",
        "tsd_len",
        "tsd_seq",
        "tsd_support",
        "n_copies_checked",
        "n_copies_with_tsd",
        "n_distinct_tsd_hits",
        "source",
        "tsd_len_support",
        "tsd_exact_support",
        "tsd_confidence",
        "n_copies_supporting_len",
        "n_copies_supporting_exact",
    ]

    with open(outfile, "w") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        for call in calls:
            row = dict(call.__dict__)

            for k, v in row.items():
                if v == "":
                    row[k] = "NA"
            
            writer.writerow(row)
