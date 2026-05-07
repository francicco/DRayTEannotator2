import argparse

from pathlib import Path
from .classify import classify_family
from .features import build_families_from_evidence
from .hmmer import parse_domtblout
from .pipeline import run_domain_annotation
from .io import load_families_tsv, write_classification_tsv, write_evidence_tsv
from .structure import load_structure_evidence_tsv
from .dfam import (
    parse_nhmmer_tblout,
    run_nhmmer,
    run_nhmmer_parallel_chunks,
)
from .mmseqs import (
    parse_mmseqs_tsv,
    run_mmseqs_rescue,
    sequence_lengths_by_clean_id,
)

from datetime import datetime
import sys


def log_step(message: str):
    now = datetime.now().isoformat(timespec="seconds")
    print(
        f"[DRayTE][classify][{now}] {message}",
        file=sys.stderr,
        flush=True,
    )

def classify_families(families):

    results = []

    for f in families:

        c = classify_family(f)

        results.append({
            "family_id": f.family_id,
            **c
        })

    return results


def main():

    parser = argparse.ArgumentParser(
        description="DRayTE evidence-based TE classifier"
    )

    input_group = parser.add_mutually_exclusive_group(required=True)

    input_group.add_argument(
        "--input",
        help="Input family evidence TSV",
    )

    input_group.add_argument(
        "--fasta",
        help="Consensus FASTA",
    )

    #
    # Protein-domain evidence: Pfam / TE-domain HMMs via hmmscan
    #

    parser.add_argument(
        "--pfam-domtblout",
        help="Existing hmmscan domtblout file with Pfam/TE protein-domain hits",
        default=None,
    )

    parser.add_argument(
        "--pfam-db",
        help="Protein HMM database for automatic hmmscan, e.g. Pfam or TE-domain HMMs",
        default=None,
    )

    parser.add_argument(
        "--pfam-outdir",
        help="Output directory for ORFs and hmmscan protein-domain results",
        default="classification_pfam",
    )

    parser.add_argument(
        "--hmmscan-bin",
        help="Path to hmmscan executable",
        default="hmmscan",
    )

    #
    # Dfam nucleotide-HMM evidence via nhmmer
    #

    parser.add_argument(
        "--dfam-db",
        help="Dfam nucleotide HMM database for nhmmer",
        default=None,
    )

    parser.add_argument(
        "--dfam-tblout",
        help="Existing nhmmer tblout file to reuse",
        default=None,
    )

    parser.add_argument(
        "--dfam-outdir",
        help="Output directory for Dfam nhmmer results",
        default="classification_dfam",
    )

    parser.add_argument(
        "--nhmmer-bin",
        help="Path to nhmmer executable",
        default="nhmmer",
    )

    parser.add_argument(
        "--chunks",
        type=int,
        default=0,
        help="Split FASTA into N chunks and run nhmmer in parallel",
    )

    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Maximum concurrent hmmscan/nhmmer jobs",
    )

    #
    # Structure evidence
    #

    parser.add_argument(
        "--structure-evidence",
        help="Structure evidence TSV",
        default=None,
    )

    #
    # MMseqs similarity rescue
    #

    parser.add_argument(
        "--mmseqs-rescue-tsv",
        help="Existing MMseqs2 easy-search TSV for similarity rescue",
        default=None,
    )

    parser.add_argument(
        "--mmseqs-bin",
        help="Path to mmseqs executable",
        default="mmseqs",
    )

    parser.add_argument(
        "--mmseqs-rescue-db",
        help="Reference FASTA for MMseqs2 rescue search",
        default=None,
    )

    parser.add_argument(
        "--mmseqs-outdir",
        help="Output directory for MMseqs2 rescue search",
        default="classification_mmseqs",
    )

    #
    # General options
    #

    parser.add_argument(
        "--output",
        required=True,
        help="Output classification TSV",
    )

    parser.add_argument(
        "--evidence-output",
        help="Optional output TSV with detailed per-family evidence",
        default=None,
    )

    parser.add_argument(
        "--min-orf-nt",
        type=int,
        default=500,
        help="Minimum ORF size in nucleotides",
    )

    parser.add_argument(
        "--forward-only",
        action="store_true",
        help="Do not scan reverse strand ORFs",
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help="Number of CPUs for hmmscan/nhmmer/MMseqs steps",
    )
    args = parser.parse_args()

    def require_db_for_missing_cache(
        evidence_name: str,
        cache_path,
        db_path,
    ):

        if db_path:
            return
    
        raise FileNotFoundError(
            f"{evidence_name} cache not found: {cache_path}. "
            f"Provide the corresponding database path to generate it automatically."
        )

    log_step("Starting classification")
    #
    # TSV mode
    #

    if args.input:
        log_step(f"Loading family evidence table: {args.input}")
        families = load_families_tsv(args.input)
        log_step(f"Loaded {len(families)} families")

    #
    # FASTA mode
    #

    else:
        log_step(f"Reading consensus FASTA: {args.fasta}")

        query_lengths = sequence_lengths_by_clean_id(
            args.fasta
        )

        log_step(f"Loaded sequence lengths for {len(query_lengths)} families")

        default_pfam_domtblout = (
            Path(args.pfam_outdir)
            / "domains.domtblout"
        )

        pfam_domtblout = (
            Path(args.pfam_domtblout)
            if args.pfam_domtblout
            else default_pfam_domtblout
        )

        if pfam_domtblout.exists():
            log_step(
                f"Parsing cached protein-domain domtblout: {pfam_domtblout}"
            )
            domain_hits = parse_domtblout(pfam_domtblout)
            log_step(f"Parsed {len(domain_hits)} protein-domain hits")

        elif args.pfam_domtblout or args.pfam_db:
            require_db_for_missing_cache(
                "Protein-domain",
                pfam_domtblout,
                args.pfam_db,
            )

            log_step(
                f"No cached protein-domain domtblout found at {pfam_domtblout}"
            )
            log_step(
                f"Running hmmscan against protein HMM database: {args.pfam_db}"
            )

            domain_hits = run_domain_annotation(
                consensus_fasta=args.fasta,
                hmm_db=args.pfam_db,
                outdir=args.pfam_outdir,
                hmmscan_bin=args.hmmscan_bin,
                min_orf_nt=args.min_orf_nt,
                include_reverse_orfs=not args.forward_only,
                cpu=args.cpu,
                chunks=args.chunks,
                max_parallel=args.jobs,
            )

            log_step(f"Detected {len(domain_hits)} protein-domain hits")

        else:
            log_step("No protein-domain evidence provided")
            domain_hits = []

        if args.structure_evidence:
            log_step(f"Loading structure evidence: {args.structure_evidence}")
            structure_evidence = load_structure_evidence_tsv(args.structure_evidence)
            log_step(f"Loaded {len(structure_evidence)} structure evidence records")
        else:
            log_step("No structure evidence provided")
            structure_evidence = []

        default_dfam_tblout = (
            Path(args.dfam_outdir)
            / "dfam.merged.tblout"
        )

        dfam_tblout = (
            Path(args.dfam_tblout)
            if args.dfam_tblout
            else default_dfam_tblout
        )

        if dfam_tblout.exists():
            log_step(f"Parsing cached Dfam tblout: {dfam_tblout}")
            dfam_hits = parse_nhmmer_tblout(dfam_tblout)
            log_step(f"Parsed {len(dfam_hits)} Dfam hits")

        elif args.dfam_tblout or args.dfam_db:
            require_db_for_missing_cache(
                "Dfam",
                dfam_tblout,
                args.dfam_db,
            )

            log_step(
                f"No cached Dfam tblout found at {dfam_tblout}"
            )
            log_step(
                f"Running nhmmer against Dfam database: {args.dfam_db}"
            )

            if args.chunks and args.chunks > 1:
                generated_tblout = run_nhmmer_parallel_chunks(
                    dfam_db=args.dfam_db,
                    consensus_fasta=args.fasta,
                    outdir=args.dfam_outdir,
                    chunks=args.chunks,
                    nhmmer_bin=args.nhmmer_bin,
                    cpu_per_job=args.cpu,
                    max_parallel=args.jobs,
                )
            else:
                generated_tblout = run_nhmmer(
                    dfam_db=args.dfam_db,
                    consensus_fasta=args.fasta,
                    tblout=dfam_tblout,
                    nhmmer_bin=args.nhmmer_bin,
                    cpu=args.cpu,
                )

            log_step(f"Parsing generated Dfam tblout: {generated_tblout}")
            dfam_hits = parse_nhmmer_tblout(generated_tblout)
            log_step(f"Parsed {len(dfam_hits)} Dfam hits")

        else:
            log_step("No Dfam evidence provided")
            dfam_hits = []

        if args.mmseqs_rescue_tsv:
            log_step(f"Parsing MMseqs rescue TSV: {args.mmseqs_rescue_tsv}")
            mmseqs_hits = parse_mmseqs_tsv(
                args.mmseqs_rescue_tsv,
                query_lengths=query_lengths,
                target_lengths=query_lengths,
                min_identity=0.85,
                min_aln_len=500,
                min_query_cov=0.40,
                min_target_cov=0.30,
            )
            log_step(f"Parsed {len(mmseqs_hits)} MMseqs rescue hits after filtering")

        elif args.mmseqs_rescue_db:
            log_step(f"Running MMseqs rescue search against {args.mmseqs_rescue_db}")
            mmseqs_tsv = run_mmseqs_rescue(
                query_fasta=args.fasta,
                reference_fasta=args.mmseqs_rescue_db,
                outdir=args.mmseqs_outdir,
                mmseqs_bin=args.mmseqs_bin,
                threads=args.cpu,
            )

            target_lengths = sequence_lengths_by_clean_id(
                args.mmseqs_rescue_db
            )

            mmseqs_hits = parse_mmseqs_tsv(
                mmseqs_tsv,
                query_lengths=query_lengths,
                target_lengths=target_lengths,
                min_identity=0.80,
                min_aln_len=500,
                min_query_cov=0.40,
                min_target_cov=0.30,
            )
            log_step(f"Parsed {len(mmseqs_hits)} MMseqs rescue hits after filtering")

        else:
            log_step("No MMseqs rescue evidence provided")
            mmseqs_hits = []

        log_step("Building per-family evidence objects")
        families = build_families_from_evidence(
            consensus_fasta=args.fasta,
            domain_hits=domain_hits,
            dfam_hits=dfam_hits,
            mmseqs_hits=mmseqs_hits,
            structure_evidence=structure_evidence,
            min_orf_nt=args.min_orf_nt,
            include_reverse_orfs=not args.forward_only,
        )
        log_step(f"Built evidence objects for {len(families)} families")

    log_step("Classifying families")
    results = classify_families(families)

    n_unknown = sum(
        1 for r in results
        if r["class"] == "Unknown"
    )

    n_classified = len(results) - n_unknown

    log_step(
        f"Classification complete: "
        f"classified={n_classified}, "
        f"unknown={n_unknown}"
    )

    log_step(
        f"Writing classification table: {args.output}"
    )

    write_classification_tsv(
        results,
        args.output
    )

    if args.evidence_output:

        log_step(
            f"Writing evidence table: "
            f"{args.evidence_output}"
        )

        write_evidence_tsv(
            families,
            results,
            args.evidence_output,
        )

    log_step(
        f"Finished: processed {len(results)} families"
    )


if __name__ == "__main__":
    main()
