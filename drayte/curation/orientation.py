from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from .repeatpeps import run_diamond_blastx, ensure_repeatpeps_db, parse_top_hit_per_query


def rewrite_repmod_header(
    input_fasta: Path,
    output_fasta: Path,
    original_name: str,
    new_header: str,
) -> Path:
    records = []
    for rec in SeqIO.parse(str(input_fasta), "fasta"):
        header = rec.id.replace(original_name, new_header)
        rec.id = header
        rec.description = header
        records.append(rec)

    with open(output_fasta, "w") as handle:
        SeqIO.write(records, handle, "fasta")
    return output_fasta


def reverse_complement_fasta(input_fasta: Path, output_fasta: Path) -> Path:
    records = []
    for rec in SeqIO.parse(str(input_fasta), "fasta"):
        rec.seq = Seq(str(rec.seq)).reverse_complement()
        records.append(rec)

    with open(output_fasta, "w") as handle:
        SeqIO.write(records, handle, "fasta")
    return output_fasta


def orient_group_files(
    family_table_tsv: Path,
    te_aid_dir: Path,
    diamond_bin: str,
    repeatpeps_db_dir: Path,
    threads: int,
    logger,
) -> None:
    db_path = ensure_repeatpeps_db(repeatpeps_db_dir, diamond_bin, logger)

    with open(family_table_tsv) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            te_class = row["top_orf_class"]
            if te_class not in {"LINE", "SINE", "LTR", "RC", "DNA"}:
                continue

            consname = row["name"]
            family = row["family"]
            group_dir = te_aid_dir / te_class

            rep_fa = group_dir / f"{consname}_rep.fa"
            msa_fa = group_dir / f"{consname}_MSA_extended.fa"
            rep_mod_fa = group_dir / f"{consname}_rep_mod.fa"
            blastx_out = group_dir / f"{consname}_extended_rep_blastx.out"

            if not rep_fa.exists():
                logger.warning("Missing rep file for orientation check: %s", rep_fa)
                continue

            clean_name = consname.replace("-rnd-", ".").replace("_family-", ".")
            new_header = f"{clean_name}#{row['class']}/{family}"

            if not rep_mod_fa.exists():
                rewrite_repmod_header(
                    input_fasta=rep_fa,
                    output_fasta=rep_mod_fa,
                    original_name=consname[:-1] if consname else consname,
                    new_header=new_header,
                )

            run_diamond_blastx(rep_fa, db_path, blastx_out, diamond_bin, threads, logger)
            top = parse_top_hit_per_query(blastx_out)

            if not top:
                logger.info("No blastx hit for %s; leaving orientation unchanged", consname)
                continue

            hit = next(iter(top.values()))
            if hit["qstart"] > hit["qend"]:
                logger.info("Reverse-complementing %s", consname)

                reverse_complement_fasta(rep_fa, rep_fa.with_suffix(".tmp"))
                rep_fa.with_suffix(".tmp").replace(rep_fa)

                if msa_fa.exists():
                    reverse_complement_fasta(msa_fa, msa_fa.with_suffix(".tmp"))
                    msa_fa.with_suffix(".tmp").replace(msa_fa)

                if rep_mod_fa.exists():
                    reverse_complement_fasta(rep_mod_fa, rep_mod_fa.with_suffix(".tmp"))
                    rep_mod_fa.with_suffix(".tmp").replace(rep_mod_fa)
