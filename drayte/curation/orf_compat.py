from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class OrfCall:
    seq_id: str
    orf_index: int
    strand: str
    frame: int
    start_nt_1based: int
    end_nt_1based: int
    nt_length: int
    aa_seq: Seq

    @property
    def fasta_id(self) -> str:
        return f"{self.seq_id}_orf{self.orf_index}"

    @property
    def fasta_description(self) -> str:
        return (
            f"strand={self.strand} frame={self.frame} "
            f"start={self.start_nt_1based} end={self.end_nt_1based} "
            f"nt_len={self.nt_length}"
        )


def _scan_stop_to_stop_orfs_single_strand(
    seq: Seq,
    seq_id: str,
    strand: str,
    minsize_nt: int,
    reverse_original_len: int | None = None,
) -> List[OrfCall]:
    """
    EMBOSS-like stop-to-stop ORF scan.
    For each frame, ORFs are segments between stop codons.
    Terminal partial ORFs are retained if they reach the sequence end and pass minsize_nt.
    Coordinates are reported in original-sequence 1-based space.
    """
    results: List[OrfCall] = []
    seq_len = len(seq)
    orf_index = 0

    for frame in range(3):
        last_stop_end = frame

        i = frame
        while i <= seq_len - 3:
            codon = str(seq[i:i+3]).upper()

            if codon in STOP_CODONS:
                orf_start = last_stop_end
                orf_end = i  # exclusive of stop codon
                nt_len = orf_end - orf_start

                if nt_len >= minsize_nt:
                    nt_seq = seq[orf_start:orf_end]
                    aa_seq = nt_seq.translate(to_stop=False)

                    if aa_seq.endswith("*"):
                        aa_seq = aa_seq[:-1]

                    orf_index += 1

                    if strand == "+":
                        start_1 = orf_start + 1
                        end_1 = orf_end
                    else:
                        assert reverse_original_len is not None
                        # map reverse-complement coordinates back to original sequence
                        start_1 = reverse_original_len - orf_end + 1
                        end_1 = reverse_original_len - orf_start

                    results.append(
                        OrfCall(
                            seq_id=seq_id,
                            orf_index=orf_index,
                            strand=strand,
                            frame=frame + 1,
                            start_nt_1based=start_1,
                            end_nt_1based=end_1,
                            nt_length=nt_len,
                            aa_seq=aa_seq,
                        )
                    )

                last_stop_end = i + 3

            i += 3

        # terminal partial ORF
        orf_start = last_stop_end
        orf_end = seq_len
        nt_len = orf_end - orf_start

        if nt_len >= minsize_nt:
            nt_seq = seq[orf_start:orf_end]
            # trim trailing partial codon
            usable = (len(nt_seq) // 3) * 3
            nt_seq = nt_seq[:usable]

            if len(nt_seq) >= minsize_nt:
                aa_seq = nt_seq.translate(to_stop=False)
                if aa_seq.endswith("*"):
                    aa_seq = aa_seq[:-1]

                orf_index += 1

                if strand == "+":
                    start_1 = orf_start + 1
                    end_1 = orf_start + len(nt_seq)
                else:
                    assert reverse_original_len is not None
                    start_1 = reverse_original_len - (orf_start + len(nt_seq)) + 1
                    end_1 = reverse_original_len - orf_start

                results.append(
                    OrfCall(
                        seq_id=seq_id,
                        orf_index=orf_index,
                        strand=strand,
                        frame=frame + 1,
                        start_nt_1based=start_1,
                        end_nt_1based=end_1,
                        nt_length=len(nt_seq),
                        aa_seq=aa_seq,
                    )
                )

    return results


def find_orfs_getorf_compatible(
    seq_record: SeqRecord,
    minsize_nt: int = 500,
    include_reverse: bool = True,
) -> List[OrfCall]:
    seq = seq_record.seq.upper()
    seq_id = seq_record.id

    results = _scan_stop_to_stop_orfs_single_strand(
        seq=seq,
        seq_id=seq_id,
        strand="+",
        minsize_nt=minsize_nt,
    )

    if include_reverse:
        rev = seq.reverse_complement()
        results.extend(
            _scan_stop_to_stop_orfs_single_strand(
                seq=rev,
                seq_id=seq_id,
                strand="-",
                minsize_nt=minsize_nt,
                reverse_original_len=len(seq),
            )
        )

    # Stable sort: longest ORFs first, then by coordinates, then strand/frame
    results.sort(
        key=lambda x: (
            -x.nt_length,
            x.start_nt_1based,
            x.end_nt_1based,
            x.strand,
            x.frame,
        )
    )

    # Renumber after sort for deterministic naming
    renumbered: List[OrfCall] = []
    for idx, call in enumerate(results, start=1):
        renumbered.append(
            OrfCall(
                seq_id=call.seq_id,
                orf_index=idx,
                strand=call.strand,
                frame=call.frame,
                start_nt_1based=call.start_nt_1based,
                end_nt_1based=call.end_nt_1based,
                nt_length=call.nt_length,
                aa_seq=call.aa_seq,
            )
        )

    return renumbered


def write_orfs_fasta(
    input_fasta: Path,
    output_fasta: Path,
    minsize_nt: int = 500,
    include_reverse: bool = True,
) -> int:
    all_records: List[SeqRecord] = []

    for rec in SeqIO.parse(str(input_fasta), "fasta"):
        calls = find_orfs_getorf_compatible(
            seq_record=rec,
            minsize_nt=minsize_nt,
            include_reverse=include_reverse,
        )

        for call in calls:
            if len(call.aa_seq) == 0:
                continue

            all_records.append(
                SeqRecord(
                    call.aa_seq,
                    id=call.fasta_id,
                    description=call.fasta_description,
                )
            )

    with open(output_fasta, "w") as handle:
        SeqIO.write(all_records, handle, "fasta")

    return len(all_records)
