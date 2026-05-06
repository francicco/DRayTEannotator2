from __future__ import annotations

from pathlib import Path
from typing import Dict

from Bio import SeqIO

from .models import Family
from .orfs import extract_orf_calls_by_family, summarize_orfs
from .domainmap import normalize_domains
from .hmmer import DomainHit, summarize_domains_by_family
from .structure import StructureEvidence, summarize_structure_evidence


def consensus_lengths(input_fasta: str | Path) -> Dict[str, int]:
    return {
        rec.id: len(rec.seq)
        for rec in SeqIO.parse(str(input_fasta), "fasta")
    }


def build_families_from_evidence(
    consensus_fasta: str | Path,
    domain_hits: list[DomainHit] | None = None,
    structure_evidence: list[StructureEvidence] | None = None,
    min_orf_nt: int = 500,
    include_reverse_orfs: bool = True,
) -> list[Family]:
    lengths = consensus_lengths(consensus_fasta)

    orf_calls = extract_orf_calls_by_family(
        consensus_fasta,
        minsize_nt=min_orf_nt,
        include_reverse=include_reverse_orfs,
    )
    orf_summary = summarize_orfs(orf_calls)

    domain_summary = summarize_domains_by_family(domain_hits or [])
    structure_summary = summarize_structure_evidence(structure_evidence or [])

    families = []

    for family_id, seq_len in lengths.items():
        raw_domains = set(domain_summary.get(family_id, {}).get("domains", []))
        domains = normalize_domains(raw_domains)
        orfs = orf_summary.get(
            family_id,
            {
                "orf_count": 0,
                "orf_max_len": 0,
            },
        )

        struct = structure_summary.get(family_id, {})

        families.append(
            Family(
                family_id=family_id,
                consensus_len=seq_len,
                n_copies=0,
                domains=domains,
                ltr_present=struct.get("ltr_present", False),
                tir_present=struct.get("tir_present", False),
                helitron_signal=struct.get("helitron_signal", False),
                tsd_present=struct.get("tsd_present", False),
                polyA_present=struct.get("polyA_present", False),
                orf_count=orfs["orf_count"],
                orf_max_len=orfs["orf_max_len"],
            )
        )

    return families
