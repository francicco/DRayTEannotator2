from __future__ import annotations

from pathlib import Path
from typing import Dict

from Bio import SeqIO

from .domainmap import normalize_domains
from .hmmer import DomainHit, summarize_domains_by_family
from .homology import homology_from_repeatmasker_header
from .ids import clean_family_id
from .models import Family
from .orfs import extract_orf_calls_by_family, summarize_orfs
from .structure import StructureEvidence, summarize_structure_evidence


def consensus_lengths(input_fasta: str | Path) -> Dict[str, int]:
    return {
        clean_family_id(rec.id): len(rec.seq)
        for rec in SeqIO.parse(str(input_fasta), "fasta")
    }


def raw_ids_by_clean_id(input_fasta: str | Path) -> Dict[str, str]:
    return {
        clean_family_id(rec.id): rec.id
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
    raw_ids = raw_ids_by_clean_id(consensus_fasta)

    orf_calls_raw = extract_orf_calls_by_family(
        consensus_fasta,
        minsize_nt=min_orf_nt,
        include_reverse=include_reverse_orfs,
    )

    orf_calls = {
        clean_family_id(k): v
        for k, v in orf_calls_raw.items()
    }

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

        hom = homology_from_repeatmasker_header(
            raw_ids.get(family_id, family_id)
        )

        if hom is None:
            homology_class = "Unknown"
            homology_superfamily = "Unknown"
            homology_score = 0.0
        else:
            homology_class = hom.homology_class
            homology_superfamily = hom.homology_superfamily
            homology_score = hom.homology_score

        families.append(
            Family(
                family_id=family_id,
                consensus_len=seq_len,
                n_copies=0,
                homology_class=homology_class,
                homology_superfamily=homology_superfamily,
                homology_score=homology_score,
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
