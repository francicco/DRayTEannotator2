import argparse

from .classify import classify_family
from .dfam import parse_nhmmer_tblout, run_nhmmer
from .features import build_families_from_evidence
from .hmmer import parse_domtblout
from .pipeline import run_domain_annotation
from .io import load_families_tsv, write_classification_tsv, write_evidence_tsv
from .structure import load_structure_evidence_tsv


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
        help="Input family TSV"
    )

    input_group.add_argument(
        "--fasta",
        help="Consensus FASTA"
    )

    parser.add_argument(
        "--domtblout",
        help="hmmscan domtblout file",
        default=None,
    )

    parser.add_argument(
        "--hmm-db",
        help="HMM database for automatic hmmscan",
        default=None,
    )

    parser.add_argument(
        "--hmmscan-bin",
        help="Path to hmmscan executable",
        default="hmmscan",
    )

    parser.add_argument(
        "--domain-outdir",
        help="Output directory for ORFs and hmmscan results",
        default="classification_domains",
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help="Number of CPUs for hmmscan",
    )

    parser.add_argument(
        "--structure-evidence",
        help="Structure evidence TSV",
        default=None,
    )

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
        "--output",
        required=True,
        help="Output classification TSV"
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
        help="Minimum ORF size in nucleotides"
    )

    parser.add_argument(
        "--forward-only",
        action="store_true",
        help="Do not scan reverse strand ORFs"
    )

    args = parser.parse_args()

    #
    # TSV mode
    #

    if args.input:

        families = load_families_tsv(args.input)

    #
    # FASTA mode
    #

    else:

        if args.domtblout:

            domain_hits = parse_domtblout(
                args.domtblout
            )

        elif args.hmm_db:

            domain_hits = run_domain_annotation(
                consensus_fasta=args.fasta,
                hmm_db=args.hmm_db,
                outdir=args.domain_outdir,
                hmmscan_bin=args.hmmscan_bin,
                min_orf_nt=args.min_orf_nt,
                include_reverse_orfs=not args.forward_only,
                cpu=args.cpu,
            )

        else:

            domain_hits = []

        if args.structure_evidence:
            structure_evidence = load_structure_evidence_tsv(
                args.structure_evidence
            )
        else:
            structure_evidence = []

        if args.dfam_tblout:
            dfam_hits = parse_nhmmer_tblout(args.dfam_tblout)

        elif args.dfam_db:
            from pathlib import Path

            dfam_tblout = (
                Path(args.dfam_outdir)
                / "dfam.tblout"
            )

            run_nhmmer(
                dfam_db=args.dfam_db,
                consensus_fasta=args.fasta,
                tblout=dfam_tblout,
                nhmmer_bin=args.nhmmer_bin,
                cpu=args.cpu,
            )

            dfam_hits = parse_nhmmer_tblout(dfam_tblout)

        else:
            dfam_hits = []

        families = build_families_from_evidence(
            consensus_fasta=args.fasta,
            domain_hits=domain_hits,
            dfam_hits=dfam_hits,
            structure_evidence=structure_evidence,
            min_orf_nt=args.min_orf_nt,
            include_reverse_orfs=not args.forward_only,
        )

    results = classify_families(families)

    write_classification_tsv(
        results,
        args.output
    )

    if args.evidence_output:
        write_evidence_tsv(
            families,
            results,
            args.evidence_output,
        )

    print(f"Classified {len(results)} families")


if __name__ == "__main__":
    main()
