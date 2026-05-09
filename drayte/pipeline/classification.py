from __future__ import annotations

import json
from pathlib import Path

from drayte.classification.features import build_families_from_evidence
from drayte.classification.run import classify_families
from drayte.classification.io import (
    write_classification_tsv,
    write_evidence_tsv,
)
from drayte.classification.structure_detect import (
    detect_tirs_from_fasta,
    write_tir_structure_tsv,
)
from drayte.classification.structure import load_structure_evidence_tsv
from drayte.classification.hmmer import parse_domtblout
from drayte.classification.pipeline import run_domain_annotation
from drayte.classification.dfam import (
    parse_nhmmer_tblout,
    run_nhmmer,
    run_nhmmer_parallel_chunks,
)
from drayte.classification.mmseqs import (
    parse_mmseqs_tsv,
    run_mmseqs_rescue,
    sequence_lengths_by_clean_id,
)
from drayte.classification.write_library import rewrite_fasta_headers
from drayte.utils.paths import ensure_dir, stage_dir
from drayte.classification.dfammap import parse_dfam_hmm_metadata
from drayte.classification.nonauto import (
    infer_nonautonomous_candidates,
)

def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def validate_file(path: Path, label: str) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")


def require_db_for_missing_cache(
    evidence_name: str,
    cache_path: Path,
    db_path: str | None,
) -> None:
    if db_path:
        return

    raise FileNotFoundError(
        f"{evidence_name} cache not found: {cache_path}. "
        f"Provide the corresponding database path to generate it automatically."
    )


def run(
    config,
    curation_result: dict,
    logger,
    stage_name: str = "classification",
    final_mode: bool = False,
    annotation_result: dict | None = None,
    refinement_result: dict | None = None,
) -> dict:
    species = config.species

    outdir = ensure_dir(stage_dir(config.outdir_path, stage_name))
    manifest = outdir / f"{stage_name}.manifest.json"

    logger.info("=" * 80)
    logger.info("STAGE: %s", stage_name)
    logger.info("Output directory: %s", outdir)
    logger.info("=" * 80)

    if manifest.exists():
        logger.info("Classification manifest already exists; skipping")
        return json.loads(manifest.read_text())

    library = Path(curation_result["final_library"])
    validate_file(library, "curated TE library")

    classification_cfg = config.extra.get("classification", {})

    min_orf_nt = int(classification_cfg.get("min_orf_nt", 300))
    forward_only = bool(classification_cfg.get("forward_only", False))

    chunks = int(classification_cfg.get("chunks", config.extra.get("chunks", 0)))
    jobs = classification_cfg.get("jobs", config.extra.get("jobs", None))
    if jobs is not None:
        jobs = int(jobs)

    cpu = int(classification_cfg.get("cpu", 1))

    #
    # Structure evidence
    #

    structure_tsv = outdir / f"{species}.structure.tsv"

    if structure_tsv.exists() and structure_tsv.stat().st_size > 0:
        logger.info("Parsing cached structure evidence: %s", structure_tsv)
        structure_evidence = load_structure_evidence_tsv(structure_tsv)
    else:
        logger.info("Detecting TIR structure evidence")
        tir_detections = detect_tirs_from_fasta(
            library,
            logger=logger,
        )
        write_tir_structure_tsv(tir_detections, structure_tsv)
        structure_evidence = load_structure_evidence_tsv(structure_tsv)

    logger.info("Loaded %d structure evidence records", len(structure_evidence))

    #
    # Pfam / protein-domain evidence
    #

    pfam_outdir = ensure_dir(outdir / "pfam")
    pfam_domtblout = pfam_outdir / "domains.domtblout"
    pfam_db = classification_cfg.get("pfam_db")

    if pfam_domtblout.exists() and pfam_domtblout.stat().st_size > 0:
        logger.info("Parsing cached protein-domain domtblout: %s", pfam_domtblout)
        domain_hits = parse_domtblout(pfam_domtblout)

    elif pfam_db:
        logger.info("Running hmmscan against protein HMM database: %s", pfam_db)
        domain_hits = run_domain_annotation(
            consensus_fasta=library,
            hmm_db=pfam_db,
            outdir=pfam_outdir,
            hmmscan_bin=classification_cfg.get("hmmscan_bin", "hmmscan"),
            min_orf_nt=min_orf_nt,
            include_reverse_orfs=not forward_only,
            cpu=cpu,
            chunks=chunks,
            max_parallel=jobs,
        )

    else:
        logger.info("No protein-domain evidence provided")
        domain_hits = []

    logger.info("Loaded %d protein-domain hits", len(domain_hits))

    #
    # Dfam evidence
    #

    dfam_outdir = ensure_dir(outdir / "dfam")
    dfam_tblout = dfam_outdir / "dfam.merged.tblout"
    dfam_db = classification_cfg.get("dfam_db")

    dfam_metadata = {}

    if dfam_db:
        logger.info("Parsing Dfam HMM metadata: %s", dfam_db)
        dfam_metadata = parse_dfam_hmm_metadata(dfam_db)
        logger.info(
            "Loaded Dfam metadata for %d models",
            len(dfam_metadata),
        )

    if dfam_tblout.exists() and dfam_tblout.stat().st_size > 0:
        logger.info("Parsing cached Dfam tblout: %s", dfam_tblout)

        dfam_hits = parse_nhmmer_tblout(
            dfam_tblout,
            dfam_metadata=dfam_metadata,
        )

    elif dfam_db:
        logger.info("Running nhmmer against Dfam database: %s", dfam_db)

        if chunks and chunks > 1:
            generated_tblout = run_nhmmer_parallel_chunks(
                dfam_db=dfam_db,
                consensus_fasta=library,
                outdir=dfam_outdir,
                chunks=chunks,
                nhmmer_bin=classification_cfg.get("nhmmer_bin", "nhmmer"),
                cpu_per_job=cpu,
                max_parallel=jobs,
            )
        else:
            generated_tblout = run_nhmmer(
                dfam_db=dfam_db,
                consensus_fasta=library,
                tblout=dfam_tblout,
                nhmmer_bin=classification_cfg.get("nhmmer_bin", "nhmmer"),
                cpu=cpu,
            )

        dfam_hits = parse_nhmmer_tblout(
            generated_tblout,
            dfam_metadata=dfam_metadata,
        )

    else:
        logger.info("No Dfam evidence provided")
        dfam_hits = []

    logger.info("Loaded %d Dfam hits", len(dfam_hits))

    #
    # MMseqs rescue evidence
    #

    mmseqs_outdir = ensure_dir(outdir / "mmseqs")
    mmseqs_tsv = mmseqs_outdir / "mmseqs_rescue.tsv"
    mmseqs_rescue_db = classification_cfg.get("mmseqs_rescue_db")

    query_lengths = sequence_lengths_by_clean_id(library)

    if mmseqs_tsv.exists() and mmseqs_tsv.stat().st_size > 0:
        logger.info("Parsing cached MMseqs rescue TSV: %s", mmseqs_tsv)

        target_lengths = sequence_lengths_by_clean_id(
            Path(mmseqs_rescue_db) if mmseqs_rescue_db else library
        )

        mmseqs_hits = parse_mmseqs_tsv(
            mmseqs_tsv,
            query_lengths=query_lengths,
            target_lengths=target_lengths,
            min_identity=float(classification_cfg.get("mmseqs_min_identity", 0.80)),
            min_aln_len=int(classification_cfg.get("mmseqs_min_aln_len", 500)),
            min_query_cov=float(classification_cfg.get("mmseqs_min_query_cov", 0.40)),
            min_target_cov=float(classification_cfg.get("mmseqs_min_target_cov", 0.30)),
        )

    elif mmseqs_rescue_db:
        logger.info("Running MMseqs rescue search against %s", mmseqs_rescue_db)

        generated_tsv = run_mmseqs_rescue(
            query_fasta=library,
            reference_fasta=mmseqs_rescue_db,
            outdir=mmseqs_outdir,
            mmseqs_bin=classification_cfg.get("mmseqs_bin", "mmseqs"),
            threads=config.threads,
        )

        target_lengths = sequence_lengths_by_clean_id(mmseqs_rescue_db)

        mmseqs_hits = parse_mmseqs_tsv(
            generated_tsv,
            query_lengths=query_lengths,
            target_lengths=target_lengths,
            min_identity=float(classification_cfg.get("mmseqs_min_identity", 0.80)),
            min_aln_len=int(classification_cfg.get("mmseqs_min_aln_len", 500)),
            min_query_cov=float(classification_cfg.get("mmseqs_min_query_cov", 0.40)),
            min_target_cov=float(classification_cfg.get("mmseqs_min_target_cov", 0.30)),
        )

    else:
        logger.info("No MMseqs rescue evidence provided")
        mmseqs_hits = []

    logger.info("Loaded %d MMseqs rescue hits", len(mmseqs_hits))

    #
    # Build evidence objects and classify
    #

    logger.info("Building per-family evidence objects")

    families = build_families_from_evidence(
        consensus_fasta=library,
        domain_hits=domain_hits,
        dfam_hits=dfam_hits,
        structure_evidence=structure_evidence,
        min_orf_nt=min_orf_nt,
        include_reverse_orfs=not forward_only,
        mmseqs_hits=mmseqs_hits,
    )

    logger.info("Built evidence objects for %d families", len(families))

    results = classify_families(families)

    #
    # Non-autonomous candidate inference
    #

    logger.info("Inferring non-autonomous TE candidates")

    nonauto_calls = infer_nonautonomous_candidates(families)

    logger.info(
        "Detected %d non-autonomous TE candidates",
        len(nonauto_calls),
    )

    classified = sum(1 for r in results if r["class"] != "Unknown")
    unknown = len(results) - classified

    logger.info(
        "Classification complete: classified=%d, unknown=%d",
        classified,
        unknown,
    )

    classifications_tsv = outdir / f"{species}.classifications.tsv"
    evidence_tsv = outdir / f"{species}.evidence.tsv"
    nonauto_tsv = outdir / f"{species}.nonauto_candidates.tsv"
    classified_library = outdir / f"{species}.classified.fa"

    write_classification_tsv(results, classifications_tsv)
    write_evidence_tsv(families, results, evidence_tsv)

    from drayte.classification.nonauto_io import (
        write_nonautonomous_tsv,
    )

    write_nonautonomous_tsv(
        nonauto_calls,
        nonauto_tsv,
    )

    rewrite_fasta_headers(
        input_fasta=library,
        classifications_tsv=evidence_tsv,
        output_fasta=classified_library,
        keep_unknown=True,
        taxon=classification_cfg.get("taxon", species),
    )

    result = {
        "stage": "classification",
        "outdir": str(outdir),
        "input_library": str(library),
        "structure_evidence": str(structure_tsv),
        "pfam_domtblout": str(pfam_domtblout),
        "dfam_tblout": str(dfam_tblout),
        "mmseqs_rescue_tsv": str(mmseqs_tsv),
        "classifications_tsv": str(classifications_tsv),
        "evidence_tsv": str(evidence_tsv),
        "nonauto_candidates_tsv": str(nonauto_tsv),
        "classified_library": str(classified_library),
        "final_library": str(classified_library),
        "n_families": len(results),
        "n_classified": classified,
        "n_unknown": unknown,
    }

    write_json(manifest, result)

    logger.info("Classification stage completed")
    logger.info("Classified library: %s", classified_library)

    return result
