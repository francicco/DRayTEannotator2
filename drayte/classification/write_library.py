from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO

from .ids import clean_family_id


UNKNOWN_VALUES = {"", "NA", "Unknown", "unknown", None}


def clean_value(value: str | None) -> str:
    if value in UNKNOWN_VALUES:
        return "Unknown"

    value = str(value).strip()

    if value in UNKNOWN_VALUES:
        return "Unknown"

    return value


def load_classifications(path: str | Path) -> dict[str, dict]:
    records = {}

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            family_id = clean_value(row.get("family_id"))

            if family_id == "Unknown":
                continue

            records[family_id] = row

    return records


def get_first_known(row: dict, keys: list[str]) -> str:
    for key in keys:
        value = clean_value(row.get(key))

        if value != "Unknown":
            return value

    return "Unknown"


def parse_label_from_fasta_id(seq_id: str) -> tuple[str, str]:
    """
    Extract family_id and RepeatMasker-style label from an existing FASTA ID.

    Examples
    --------
    fam1#LINE/L2 -> fam1, LINE/L2
    fam1#Unknown -> fam1, Unknown
    fam1 -> fam1, Unknown
    """
    seq_id = seq_id.split()[0]

    if "#" not in seq_id:
        return clean_family_id(seq_id), "Unknown"

    family_id, label = seq_id.split("#", 1)
    family_id = clean_family_id(family_id)
    label = clean_value(label)

    return family_id, label


def split_rm_label(label: str) -> tuple[str, str]:
    label = clean_value(label)

    if label == "Unknown":
        return "Unknown", "Unknown"

    if "/" in label:
        rm_class, superfamily = label.split("/", 1)
        return clean_value(rm_class), clean_value(superfamily)

    return clean_value(label), "Unknown"


def infer_rm_class_from_row(row: dict) -> str:
    return get_first_known(
        row,
        [
            "order",
            "final_order",
        ],
    )

def infer_superfamily_from_row(row: dict) -> str:
    return get_first_known(
        row,
        [
            "superfamily",
            "final_superfamily",
        ],
    )

def infer_final_class(row: dict) -> str:
    return get_first_known(
        row,
        [
            "class",
            "final_class",
        ],
    )

def normalize_rm_label(rm_class: str, superfamily: str) -> tuple[str, str]:
    rm_class = clean_value(rm_class)
    superfamily = clean_value(superfamily)

    if rm_class == "TIR":
        rm_class = "DNA"

    if rm_class == "Helitron":
        rm_class = "RC"
        if superfamily == "Unknown":
            superfamily = "Helitron"

    if rm_class == "Penelope":
        rm_class = "LINE"
        superfamily = "Penelope"

    sf_map = {
        "Tc1-Mariner_like": "TcMar-Mariner",
        "TcMar-Mariner": "TcMar-Mariner",
        "PIF-Harbinger_or_CACTA_like": "PIF-Harbinger",
        "CACTA_like": "PIF-Harbinger",
        "2bp_TSD_TIR_ambiguous": "Unknown",
        "BovB": "RTE-BovB",
        "CR1": "CR1",
        "L2": "L2",
        "I": "I",
    }

    if superfamily in sf_map:
        superfamily = sf_map[superfamily]

    return rm_class, superfamily

def repeatmasker_label(
    row: dict | None,
    original_label: str = "Unknown",
) -> str:
    if row is None:
        return clean_value(original_label)

    final_class = infer_final_class(row)
    rm_class = infer_rm_class_from_row(row)
    superfamily = infer_superfamily_from_row(row)
    rm_class, superfamily = normalize_rm_label(rm_class, superfamily)

    if rm_class == "Helitron":
        rm_class = "RC"
        if superfamily == "Unknown":
            superfamily = "Helitron"

    if final_class == "Unknown" or rm_class == "Unknown":
        return "Unknown"
    
    if superfamily == "Unknown":
        return rm_class
    
    return f"{rm_class}/{superfamily}"

def format_header(
    family_id: str,
    label: str,
    taxon: str | None = None,
) -> str:
    header = f"{family_id}#{label}"

    if taxon:
        header = f"{header} @{taxon}"

    return header


def rewrite_fasta_headers(
    input_fasta: str | Path,
    classifications_tsv: str | Path,
    output_fasta: str | Path,
    keep_unknown: bool = True,
    taxon: str | None = None,
) -> int:
    classifications = load_classifications(classifications_tsv)

    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0

    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(str(input_fasta), "fasta"):
            family_id, original_label = parse_label_from_fasta_id(rec.id)

            row = classifications.get(family_id)

            label = repeatmasker_label(
                row,
                original_label=original_label,
            )

            if label == "Unknown" and not keep_unknown:
                continue

            new_header = format_header(
                family_id=family_id,
                label=label,
                taxon=taxon,
            )

            rec.id = new_header
            rec.name = new_header
            rec.description = ""

            SeqIO.write(rec, out, "fasta")
            n_written += 1

    return n_written


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Rewrite TE FASTA headers using DRayTE classification results "
            "to generate a RepeatMasker-compatible classified library."
        )
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Input consensus FASTA",
    )

    parser.add_argument(
        "--classifications",
        required=True,
        help=(
            "DRayTE classification TSV or evidence TSV. "
            "Both class/order/superfamily and final_class/final_order/"
            "final_superfamily formats are supported."
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output RepeatMasker-compatible FASTA",
    )

    parser.add_argument(
        "--taxon",
        default=None,
        help="Optional taxon label appended to FASTA headers as '@taxon'",
    )

    parser.add_argument(
        "--drop-unknown",
        action="store_true",
        help="Exclude families classified as Unknown",
    )

    args = parser.parse_args()

    n_written = rewrite_fasta_headers(
        input_fasta=structure_library,
        classifications_tsv=args.classifications,
        output_fasta=args.output,
        keep_unknown=not args.drop_unknown,
        taxon=args.taxon,
    )

    print(f"Wrote {n_written} sequences to {args.output}")


if __name__ == "__main__":
    main()
