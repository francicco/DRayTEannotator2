from __future__ import annotations

import csv
import logging
from Bio.Seq import Seq
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta


LOGGER = logging.getLogger(__name__)


@dataclass
class BlastHit:
    query: str
    subject: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float

    @property
    def strand(self) -> str:
        return "+" if self.sstart <= self.send else "-"

    @property
    def start0(self) -> int:
        return min(self.sstart, self.send) - 1

    @property
    def end0(self) -> int:
        return max(self.sstart, self.send)


def read_library_ids(library_fasta: str | Path) -> List[str]:
    return [rec.id for rec in SeqIO.parse(str(library_fasta), "fasta")]


def parse_blast_outfmt6(blast_file: str | Path) -> List[BlastHit]:
    hits: List[BlastHit] = []
    with open(blast_file) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            hits.append(
                BlastHit(
                    query=row[0],
                    subject=row[1],
                    pident=float(row[2]),
                    length=int(row[3]),
                    mismatch=int(row[4]),
                    gapopen=int(row[5]),
                    qstart=int(row[6]),
                    qend=int(row[7]),
                    sstart=int(row[8]),
                    send=int(row[9]),
                    evalue=float(row[10]),
                    bitscore=float(row[11]),
                )
            )
    return hits


def group_top_hits_by_query(
    hits: List[BlastHit],
    max_hits_per_query: int = 50,
) -> Dict[str, List[BlastHit]]:
    grouped: Dict[str, List[BlastHit]] = {}
    for hit in hits:
        grouped.setdefault(hit.query, []).append(hit)

    selected: Dict[str, List[BlastHit]] = {}
    for query, qhits in grouped.items():
        qhits_sorted = sorted(
            qhits,
            key=lambda h: (-h.bitscore, h.evalue, -h.length),
        )
        selected[query] = qhits_sorted[:max_hits_per_query]

    return selected


def extract_hit_sequences(
    genome_fasta: str | Path,
    grouped_hits: Dict[str, List[BlastHit]],
    output_dir: str | Path,
    flank_left: int = 100,
    flank_right: int = 100,
) -> Dict[str, Path]:
    genome = Fasta(str(genome_fasta), as_raw=True, sequence_always_upper=True)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    outfiles: Dict[str, Path] = {}

    for query, hits in grouped_hits.items():
        outfile = output_dir / f"{query}.fa"
        records: List[SeqRecord] = []

        for idx, hit in enumerate(hits, start=1):
            chrom = hit.subject
            contig_len = len(genome[chrom])

            start = max(0, hit.start0 - flank_left)
            end = min(contig_len, hit.end0 + flank_right)

            seq = genome[chrom][start:end]
            if hit.strand == "-":
                from Bio.Seq import Seq
                seq = str(Seq(seq).reverse_complement())

            record = SeqRecord(
                seq=seq,
                id=f"{query}_hit{idx}",
                description=f"{chrom}:{start+1}-{end}({hit.strand}) bitscore={hit.bitscore}",
            )
            records.append(record)

        with open(outfile, "w") as handle:
            SeqIO.write(records, handle, "fasta")

        outfiles[query] = outfile
        LOGGER.info("Wrote %d extracted hits for %s -> %s", len(records), query, outfile)

    return outfiles


def run_extract_align(
    genome_fasta: str | Path,
    blast_file: str | Path,
    library_fasta: str | Path,
    output_dir: str | Path,
    max_hits_per_query: int = 50,
    flank_left: int = 100,
    flank_right: int = 100,
) -> Dict[str, Path]:
    output_dir = Path(output_dir)
    cat_dir = output_dir / "catTEfiles"
    cat_dir.mkdir(parents=True, exist_ok=True)

    library_ids = set(read_library_ids(library_fasta))
    hits = parse_blast_outfmt6(blast_file)
    hits = [h for h in hits if h.query in library_ids]

    grouped = group_top_hits_by_query(hits, max_hits_per_query=max_hits_per_query)
    return extract_hit_sequences(
        genome_fasta=genome_fasta,
        grouped_hits=grouped,
        output_dir=cat_dir,
        flank_left=flank_left,
        flank_right=flank_right,
    )
