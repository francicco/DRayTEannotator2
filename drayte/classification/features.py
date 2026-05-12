from __future__ import annotations

from pathlib import Path
from typing import Dict

from Bio import SeqIO

from .dfam import DfamHit, best_dfam_hits_by_family
from .dfammap import infer_te_from_dfam_model
from .domainmap import normalize_domains
from .hmmer import DomainHit
from .homology import homology_from_repeatmasker_header
from .ids import clean_family_id
from .mmseqs import MMseqsHit, infer_rescue_evidence
from .models import Family
from .orfs import extract_orf_calls_by_family
from .parsers.repeatmasker import parse_repeatmasker_header
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


def header_labels_by_clean_id(input_fasta: str | Path) -> Dict[str, dict]:
    labels = {}

    for rec in SeqIO.parse(str(input_fasta), "fasta"):
        family_id = clean_family_id(rec.id)
        parsed = parse_repeatmasker_header(rec.id)

        if parsed is None:
            labels[family_id] = {
                "original_header": rec.id,
                "header_class": "Unknown",
                "header_superfamily": "Unknown",
            }
        else:
            labels[family_id] = {
                "original_header": rec.id,
                "header_class": parsed.rm_class,
                "header_superfamily": parsed.rm_superfamily,
            }

    return labels


def summarize_orfs(orf_calls_by_family):
    summary = {}

    for family_id, calls in orf_calls_by_family.items():
        if not calls:
            summary[family_id] = {
                "orf_count": 0,
                "orf_max_len": 0,
            }
            continue

        summary[family_id] = {
            "orf_count": len(calls),
            "orf_max_len": max(c.nt_length for c in calls),
        }

    return summary


def summarize_domains_by_family(domain_hits: list[DomainHit]):
    summary = {}

    for hit in domain_hits:
        fam = clean_family_id(hit.family_id)

        if fam not in summary:
            summary[fam] = {
                "domains": set(),
            }

        summary[fam]["domains"].add(hit.domain)

    return summary


def build_families_from_evidence(
    consensus_fasta: str | Path,
    domain_hits: list[DomainHit] | None = None,
    dfam_hits: list[DfamHit] | None = None,
    mmseqs_hits: list[MMseqsHit] | None = None,
    structure_evidence: list[StructureEvidence] | None = None,
    min_orf_nt: int = 500,
    include_reverse_orfs: bool = True,
) -> list[Family]:
    lengths = consensus_lengths(consensus_fasta)
    raw_ids = raw_ids_by_clean_id(consensus_fasta)
    header_labels = header_labels_by_clean_id(consensus_fasta)

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
    best_dfam = best_dfam_hits_by_family(dfam_hits or [])
    rescue_summary = infer_rescue_evidence(mmseqs_hits or [])
    structure_summary = summarize_structure_evidence(structure_evidence or [])

    tsd_summary = {}

    families = []

    # Build Family objects from the union of all evidence sources, not only
    # from sequences present in the current consensus FASTA. This is important
    # for structure-derived evidence produced from the structure input FASTA,
    # where extension-derived consensus IDs may not all be represented in the
    # classification FASTA used here.
    all_family_ids = set(lengths)
    all_family_ids.update(structure_summary)
    all_family_ids.update(domain_summary)
    all_family_ids.update(best_dfam)
    all_family_ids.update(rescue_summary)
    all_family_ids.update(orf_summary)

    for family_id in sorted(all_family_ids):
        seq_len = lengths.get(family_id)

        if seq_len is None:
            struct = structure_summary.get(family_id, {})
            seq_len = struct.get("consensus_len")

        if seq_len is None:
            seq_len = 0

        raw_domains = set(
            domain_summary.get(family_id, {}).get("domains", [])
        )
        domains = normalize_domains(raw_domains)

        orfs = orf_summary.get(
            family_id,
            {
                "orf_count": 0,
                "orf_max_len": 0,
            },
        )

        struct = structure_summary.get(family_id, {})

        header = header_labels.get(
            family_id,
            {
                "original_header": family_id,
                "header_class": "Unknown",
                "header_superfamily": "Unknown",
            },
        )

        dfam_hit = best_dfam.get(family_id)

        if dfam_hit is None:
            dfam_class = None
            dfam_order = None
            dfam_superfamily = None
            dfam_model = None
            dfam_score = 0.0
        else:
            dfam_class = getattr(dfam_hit, "dfam_class", "Unknown")
            dfam_order = getattr(dfam_hit, "dfam_order", "Unknown")
            dfam_superfamily = getattr(
                dfam_hit,
                "dfam_superfamily",
                "Unknown",
            )
            
            if (
                dfam_class in {None, "", "Unknown"}
                and dfam_order in {None, "", "Unknown"}
                and dfam_superfamily in {None, "", "Unknown"}
            ):
                (
                    dfam_class,
                    dfam_order,
                    dfam_superfamily,
                ) = infer_te_from_dfam_model(dfam_hit.model_name)

            dfam_model = dfam_hit.model_name
            dfam_score = dfam_hit.score

        rescue = rescue_summary.get(family_id)

        if rescue is None:
            rescue_class = "Unknown"
            rescue_order = "Unknown"
            rescue_superfamily = "Unknown"
            rescue_identity = 0.0
            rescue_aln_len = 0
            rescue_bits = 0.0
            rescue_target = ""
        else:
            rescue_class = rescue.rescue_class
            rescue_order = rescue.rescue_order
            rescue_superfamily = rescue.rescue_superfamily
            rescue_identity = rescue.rescue_identity
            rescue_aln_len = rescue.rescue_aln_len
            rescue_bits = rescue.rescue_bits
            rescue_target = rescue.rescue_target

        hom = homology_from_repeatmasker_header(
            raw_ids.get(family_id, family_id)
        )

        if hom is None:
            homology_class = "Unknown"
            homology_order = "Unknown"
            homology_superfamily = "Unknown"
            homology_score = 0.0
        else:
            homology_class = hom.homology_class
            homology_order = hom.homology_order
            homology_superfamily = hom.homology_superfamily
            homology_score = hom.homology_score

        families.append(
            Family(
                family_id=family_id,
                consensus_len=seq_len,
                n_copies=0,
                original_header=header["original_header"],
                header_class=header["header_class"],
                header_superfamily=header["header_superfamily"],
                homology_class=homology_class,
                homology_order=homology_order,
                homology_superfamily=homology_superfamily,
                homology_score=homology_score,
                dfam_class=dfam_class,
                dfam_order=dfam_order,
                dfam_superfamily=dfam_superfamily,
                dfam_model=dfam_model,
                dfam_score=dfam_score,
                rescue_class=rescue_class,
                rescue_order=rescue_order,
                rescue_superfamily=rescue_superfamily,
                rescue_identity=rescue_identity,
                rescue_aln_len=rescue_aln_len,
                rescue_bits=rescue_bits,
                rescue_target=rescue_target,
                domains=domains,
                tir_present=struct.get("tir_present", False),
                tir_len=struct.get("tir_len", 0),
                tir_identity=struct.get("tir_identity", 0.0),
                tir_confidence=struct.get("tir_confidence", "LOW"),
                tir_5p_motif=struct.get("tir_5p_motif", "NA"),
                tir_3p_motif=struct.get("tir_3p_motif", "NA"),
                tir_terminal_motif=struct.get("tir_terminal_motif", "NA"),
                tir_grade=struct.get("tir_grade", "ABSENT"),
                structural_type=struct.get("structural_type", "unknown"),
                tir_tsd_element=struct.get("tir_tsd_element", struct.get("mite_like_structure", False)),
                structural_superfamily_hint=struct.get("superfamily_hint", "unknown"),
                structural_superfamily_confidence=struct.get("superfamily_confidence", "NONE"),
                helitron_signal=struct.get("helitron_signal", False),
                tsd_present=struct.get("tsd_present", False),
                polyA_present=struct.get("polyA_present", False),
                polyA_tail_len=struct.get("polyA_tail_len", 0),
                polyA_fraction=struct.get("polyA_fraction", 0.0),
                poliii_a_box=struct.get("poliii_a_box", False),
                poliii_b_box=struct.get("poliii_b_box", False),
                poliii_score=struct.get("poliii_score", 0.0),
                sine_score=struct.get("sine_score", 0.0),
                sine_candidate=struct.get("sine_candidate", False),
                sine_confidence=struct.get("sine_confidence", "LOW"),
                ltr_present=struct.get("ltr_present", False),
                ltr_structural_type=struct.get("ltr_structural_type", "none"),
                ltr_score=struct.get("ltr_score", 0.0),
                ltr_confidence=struct.get("ltr_confidence", "LOW"),
                tg_ca_motif=struct.get("tg_ca_motif", False),
                ppt_like=struct.get("ppt_like", False),
                orf_count=orfs["orf_count"],
                orf_max_len=orfs["orf_max_len"],
                tsd_seq=struct.get("tsd_seq", ""),
                tsd_len=struct.get("tsd_len", 0),
                tsd_support=struct.get("tsd_support", 0.0),
            )
        )

    return families

