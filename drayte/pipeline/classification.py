from __future__ import annotations

import hashlib
import json
from pathlib import Path

from Bio import SeqIO

from drayte.classification.complexity_detect import detect_complexity_from_fasta
from drayte.classification.features import (
    build_families_from_evidence,
    consensus_lengths,
)
from drayte.classification.run import classify_families
from drayte.classification.io import (
    write_classification_tsv,
    write_evidence_tsv,
)
from drayte.classification.structure_detect import (
    detect_tirs_from_fasta,
)
from drayte.classification.structure import load_structure_evidence_tsv
from drayte.classification.structure_merge import (
    merge_structure_evidence,
    write_structure_summary,
)
from drayte.classification.tsd_from_repam import (
    infer_tsds_from_extensionwork,
    write_repam_tsd_tsv,
)
from drayte.classification.sine_detect import detect_sines_from_fasta
from drayte.classification.hmmer import parse_domtblout
from drayte.classification.pipeline import run_domain_annotation
from drayte.classification.dfam import (
    cached_nhmmer_chunk_count,
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
from drayte.classification.structure_inputs import resolve_structure_input_fasta
from drayte.utils.paths import ensure_dir, stage_dir
from drayte.classification.dfammap import parse_dfam_hmm_metadata
from drayte.classification.nonauto import (
    infer_nonautonomous_candidates,
)
from drayte.classification.ltr_detect import detect_ltrs_from_fasta


CLASSIFICATION_SIGNATURE_SCHEMA_VERSION = 2
CLASSIFICATION_EXECUTION_ONLY_KEYS = frozenset({"chunks", "jobs", "cpu"})
CLASSIFICATION_CORE_SIGNATURE_KEYS = (
    "stage_name",
    "final_mode",
    "library_sha256",
    "structure_library_sha256",
    "extensionwork_sha256",
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def selected_extensionwork_sha256(extensionwork_dir: Path) -> str:
    digest = hashlib.sha256()
    if not extensionwork_dir.exists():
        digest.update(b"missing")
        return digest.hexdigest()

    selected: list[Path] = []
    for filename in ("repam-ranges.tsv", "repam-repseq.fa"):
        selected.extend(extensionwork_dir.glob(f"*/{filename}"))

    for path in sorted(
        selected,
        key=lambda item: str(item.relative_to(extensionwork_dir)),
    ):
        relative = str(path.relative_to(extensionwork_dir))
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)

    return digest.hexdigest()


def classification_input_signature(
    library: Path,
    structure_library: Path,
    extensionwork_dir: Path,
    classification_cfg: dict,
    stage_name: str,
    final_mode: bool,
) -> dict:
    result_config = {
        key: value
        for key, value in classification_cfg.items()
        if key not in CLASSIFICATION_EXECUTION_ONLY_KEYS
    }
    execution_config = {
        key: classification_cfg[key]
        for key in sorted(CLASSIFICATION_EXECUTION_ONLY_KEYS)
        if key in classification_cfg
    }
    result_config_json = json.dumps(result_config, sort_keys=True, default=str)
    return {
        "schema_version": CLASSIFICATION_SIGNATURE_SCHEMA_VERSION,
        "stage_name": stage_name,
        "final_mode": final_mode,
        "library_sha256": sha256_file(library),
        "structure_library_sha256": sha256_file(structure_library),
        "extensionwork_sha256": selected_extensionwork_sha256(extensionwork_dir),
        "result_config_sha256": hashlib.sha256(
            result_config_json.encode("utf-8")
        ).hexdigest(),
        "execution_config": execution_config,
    }


def classification_signatures_compatible(
    cached: dict,
    current: dict,
) -> tuple[bool, bool]:
    """
    Compare classification checkpoints without treating scheduler settings as
    biological inputs.

    Returns ``(compatible, legacy_migration)``.  Version-one signatures hashed
    all configuration fields together, so they cannot distinguish a scheduling
    change from a result-affecting change.  They may be migrated only when all
    biological input digests still match; the caller additionally validates the
    classified library and every required output before reuse.
    """
    if any(cached.get(key) != current.get(key) for key in CLASSIFICATION_CORE_SIGNATURE_KEYS):
        return False, False

    cached_version = cached.get("schema_version", 1)
    if cached_version == 1:
        return True, True

    if cached_version != CLASSIFICATION_SIGNATURE_SCHEMA_VERSION:
        return False, False

    return (
        cached.get("result_config_sha256") == current.get("result_config_sha256"),
        False,
    )


def manifest_outputs_are_valid(result: dict) -> bool:
    required_keys = (
        "structure_evidence",
        "classifications_tsv",
        "evidence_tsv",
        "classified_library",
    )
    for key in required_keys:
        value = result.get(key)
        if not value:
            return False
        path = Path(value)
        if not path.exists() or path.stat().st_size == 0:
            return False
    return True


def fasta_sequences_by_family(fasta: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for record in SeqIO.parse(str(fasta), "fasta"):
        family_id = record.id.split("#", 1)[0]
        if family_id in sequences:
            raise ValueError(f"Duplicate family identifier in {fasta}: {family_id}")
        sequences[family_id] = str(record.seq).upper()
    return sequences


def classified_library_matches_input(input_fasta: Path, output_fasta: Path) -> bool:
    if not output_fasta.exists() or output_fasta.stat().st_size == 0:
        return False
    try:
        input_sequences = fasta_sequences_by_family(input_fasta)
        output_sequences = fasta_sequences_by_family(output_fasta)
    except (OSError, ValueError):
        return False
    return bool(input_sequences) and input_sequences == output_sequences


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

    library = Path(curation_result["final_library"])
    validate_file(library, "curated TE library")

    classification_cfg = config.extra.get("classification", {})

    extensionwork_dir = Path(
        classification_cfg.get(
            "repam_extensionwork_dir",
            config.outdir_path / "extension" / "extensionwork",
        )
    )

    # Structural evidence is generated for the broad discovery/extension
    # family universe, not only the curated final library. The curated library
    # can omit families that still have REPAM/TSD evidence in extensionwork.
    structure_library = resolve_structure_input_fasta(
        config=config,
        classification_cfg=classification_cfg,
        outdir=outdir,
        fallback_fasta=library,
        logger=logger,
    )
    validate_file(structure_library, "structure input library")

    input_signature = classification_input_signature(
        library=library,
        structure_library=structure_library,
        extensionwork_dir=extensionwork_dir,
        classification_cfg=classification_cfg,
        stage_name=stage_name,
        final_mode=final_mode,
    )

    if manifest.exists():
        try:
            cached_result = json.loads(manifest.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning("Could not read classification manifest; resuming: %s", exc)
        else:
            cached_signature = cached_result.get("input_signature")
            legacy_signature_migration = False
            if cached_signature is not None:
                signatures_compatible, legacy_signature_migration = (
                    classification_signatures_compatible(
                        cached_signature,
                        input_signature,
                    )
                )
                if not signatures_compatible:
                    raise RuntimeError(
                        "Existing classification outputs were created from different "
                        "biological inputs or result-affecting settings. Move the "
                        "existing classification stage directory aside, or use a new "
                        "outdir, before restarting."
                    )

            if cached_signature is None:
                output_value = cached_result.get("classified_library")
                if not output_value or not classified_library_matches_input(
                    structure_library,
                    Path(output_value),
                ):
                    raise RuntimeError(
                        "Legacy classification manifest cannot be matched to the "
                        "current input library. Move the existing classification "
                        "stage directory aside, or use a new outdir, before "
                        "restarting."
                    )

            if legacy_signature_migration:
                output_value = cached_result.get("classified_library")
                if not output_value or not classified_library_matches_input(
                    structure_library,
                    Path(output_value),
                ):
                    raise RuntimeError(
                        "Legacy classification checkpoint cannot be safely migrated: "
                        "the classified library does not match the current structure "
                        "input library. Move the existing classification stage "
                        "directory aside, or use a new outdir, before restarting."
                    )

            if manifest_outputs_are_valid(cached_result):
                if legacy_signature_migration:
                    logger.warning(
                        "Migrating legacy classification checkpoint after validating "
                        "unchanged biological inputs and outputs; execution-only "
                        "settings such as chunks, jobs, and cpu no longer invalidate "
                        "completed classification results"
                    )
                cached_result["input_signature"] = input_signature
                write_json(manifest, cached_result)
                logger.info("Validated existing classification outputs; skipping")
                return cached_result

            logger.warning("Classification cache is incomplete; resuming stage")

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
    # Build one merged structure_summary file that contains both terminal TIR
    # calls and REPAM-derived TSD support. load_structure_evidence_tsv() can
    # read this wide table and convert it into evidence objects for the
    # downstream Family model.

    structure_tsv = outdir / f"{species}.structure_summary.tsv"
    repam_tsd_tsv = outdir / f"{species}.repam_tsd.tsv"

    legacy_structure_tsv = outdir / f"{species}.structure.tsv"

    logger.info("Detecting low-complexity/simple-repeat evidence")
    complexity_detections = detect_complexity_from_fasta(structure_library)

    if structure_tsv.exists() and structure_tsv.stat().st_size > 0:
        logger.info("Parsing cached merged structure evidence: %s", structure_tsv)
        structure_evidence = load_structure_evidence_tsv(structure_tsv)

    elif legacy_structure_tsv.exists() and legacy_structure_tsv.stat().st_size > 0:
        logger.info(
            "Parsing legacy structure evidence: %s",
            legacy_structure_tsv,
        )
        structure_evidence = load_structure_evidence_tsv(legacy_structure_tsv)
        structure_tsv = legacy_structure_tsv

    else:
        logger.info("Detecting terminal TIR structure evidence")
        tir_detections = detect_tirs_from_fasta(
            structure_library,
            threads=config.threads,
            logger=logger,
        )

        logger.info("Detecting first-pass structural SINE evidence")
        sine_detections = detect_sines_from_fasta(
            structure_library,
            max_len=int(classification_cfg.get("sine_max_len", 700)),
            tail_window=int(classification_cfg.get("sine_tail_window", 80)),
            head_window=int(classification_cfg.get("sine_head_window", 120)),
            min_poly_run=int(classification_cfg.get("sine_min_poly_run", 10)),
            min_tail_at_fraction=float(classification_cfg.get("sine_min_tail_at_fraction", 0.65)),
            min_score=float(classification_cfg.get("sine_min_score", 0.70)),
        )

        logger.info("Detecting terminal LTR/TRIM/LARD structure evidence")
        ltr_detections = detect_ltrs_from_fasta(structure_library)

        if extensionwork_dir.exists():
            logger.info("Inferring REPAM TSD evidence from %s", extensionwork_dir)
            tsd_calls = infer_tsds_from_extensionwork(
                extensionwork_dir=extensionwork_dir,
                min_len=int(classification_cfg.get("tsd_min_len", 2)),
                max_len=int(classification_cfg.get("tsd_max_len", 15)),
                min_copies=int(classification_cfg.get("tsd_min_copies", 3)),
                min_support=float(classification_cfg.get("tsd_min_support", 0.30)),
            )
            write_repam_tsd_tsv(tsd_calls, repam_tsd_tsv)
        else:
            logger.info(
                "No REPAM extensionwork directory found; TSD evidence unavailable: %s",
                extensionwork_dir,
            )
            tsd_calls = []

        merged_structure = merge_structure_evidence(
            complexity_calls=complexity_detections,
            tir_calls=tir_detections,
            tsd_calls=tsd_calls,
            sine_calls=sine_detections,
            ltr_calls=ltr_detections,
            consensus_lengths=consensus_lengths(structure_library),
        )
        write_structure_summary(merged_structure, structure_tsv)
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

        cached_dfam_chunks = cached_nhmmer_chunk_count(dfam_outdir)
        if cached_dfam_chunks is not None:
            logger.info(
                "Found existing Dfam checkpoint containing %d chunks",
                cached_dfam_chunks,
            )

        if (chunks and chunks > 1) or cached_dfam_chunks is not None:
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
        consensus_fasta=structure_library,
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
        input_fasta=structure_library,
        classifications_tsv=classifications_tsv,
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
        "input_signature": input_signature,
    }

    write_json(manifest, result)

    logger.info("Classification stage completed")
    logger.info("Classified library: %s", classified_library)

    return result
