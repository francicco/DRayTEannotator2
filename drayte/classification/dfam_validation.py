"""Optional final Dfam validation using the existing NHMMER and classifier code.

Legacy flat dfam_db runs retain their historical behavior. Explicit dfam.enabled
uses this conservative policy and never appends reference sequences to a library.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
from dataclasses import asdict, replace
from pathlib import Path

from Bio import SeqIO

from .classify import classify_family, is_known
from .dfam import DfamHit, parse_nhmmer_tblout, run_nhmmer_parallel_chunks
from .dfammap import parse_dfam_hmm_metadata

POLICY_VERSION = 1
DEFAULTS = {
    "enabled": False, "curated_only": True, "engine": "nhmmer",
    "taxon": None, "ancestors": True, "descendants": True,
    "min_query_coverage": 0.10, "min_model_coverage": 0.50,
    "max_evalue": 1e-5, "min_score": 20.0, "require_ga": True,
    "competition_score_fraction": 0.95, "on_missing": "error",
    # Header and nucleotide-reference homology can already derive from Dfam.
    # An explicit provenance assertion is required for confidence promotion.
    "independent": False,
}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def options(config):
    if not isinstance(config, dict):
        raise ValueError("classification.dfam must be a mapping")
    cfg = {**DEFAULTS, **config}
    for key in ("enabled", "curated_only", "ancestors", "descendants", "require_ga", "independent"):
        if not isinstance(cfg[key], bool):
            raise ValueError(f"classification.dfam.{key} must be boolean")
    if cfg["engine"] != "nhmmer":
        raise ValueError("classification.dfam.engine must be nhmmer (nucleotide HMMs)")
    if cfg["on_missing"] not in {"error", "warn"}:
        raise ValueError("classification.dfam.on_missing must be error or warn")
    for key in ("min_query_coverage", "min_model_coverage", "competition_score_fraction"):
        value = float(cfg[key])
        if not math.isfinite(value) or not 0 < value <= 1:
            raise ValueError(f"classification.dfam.{key} must be greater than 0 and at most 1")
        cfg[key] = value
    for key in ("max_evalue", "min_score"):
        cfg[key] = float(cfg[key])
        if not math.isfinite(cfg[key]) or cfg[key] <= 0:
            raise ValueError(f"classification.dfam.{key} must be positive and finite")
    for key in ("threads", "chunks", "jobs"):
        if key in cfg and (not isinstance(cfg[key], int) or isinstance(cfg[key], bool) or cfg[key] < 1):
            raise ValueError(f"classification.dfam.{key} must be a positive integer")
    return cfg


def executable(command):
    resolved = shutil.which(str(command))
    if not resolved:
        raise FileNotFoundError(f"Dfam dependency unavailable: {command}")
    return resolved


def capture(command):
    result = subprocess.run(command, text=True, capture_output=True)
    # FamDB 3.0 prints argparse help and exits 1 rather than 0.
    help_exit = "--help" in command and result.returncode == 1 and "usage:" in result.stdout
    if result.returncode and not help_exit:
        raise subprocess.CalledProcessError(result.returncode, command, result.stdout, result.stderr)
    return result.stdout + result.stderr


def write_json(path, value):
    partial = Path(str(path) + ".partial")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(partial, path)


def prepare_reference(cfg, outdir):
    """Use the installed FamDB CLI; HDF5 layout is owned by FamDB, not us."""
    if not cfg.get("database"):
        raise FileNotFoundError("classification.dfam.database is required when enabled")
    source = Path(cfg["database"]).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"Dfam reference does not exist: {source}")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    provenance = {"source": str(source), "policy_version": POLICY_VERSION}
    if source.is_dir():
        if not cfg.get("taxon"):
            raise ValueError("A taxon is required to export a FamDB directory")
        famdb = executable(cfg.get("famdb_bin", "famdb.py"))
        help_text = capture([famdb, "families", "--help"])
        for flag in ("--format", "--curated", "--ancestors", "--descendants"):
            if flag not in help_text:
                raise ValueError(f"Installed FamDB lacks required option {flag}")
        version_help = capture([famdb, "--help"])
        info = capture([famdb, "-i", str(source), "info"])
        provenance.update(famdb_help=help_text, famdb_version_help=version_help, database_info=info,
                          famdb_sha256=sha256(famdb))
        # Only HDF5 files participate; export bytes below provide the portable
        # reference digest. Avoid rehashing multi-gigabyte partitions every run.
        partitions = [{"name": p.name, "size": p.stat().st_size,
                       "mtime_ns": p.stat().st_mtime_ns} for p in sorted(source.glob("*.h5"))]
        if not partitions:
            raise FileNotFoundError(f"No FamDB HDF5 partitions found in {source}")
        cmd = [famdb, "-i", str(source), "families", "--format", "hmm"]
        if cfg["ancestors"]:
            cmd.append("--ancestors")
        if cfg["descendants"]:
            cmd.append("--descendants")
        if cfg["curated_only"]:
            cmd.append("--curated")
        cmd.append(str(cfg["taxon"]))
        signature = {"partitions": partitions, "command": cmd, "info": info,
                     "famdb_sha256": provenance["famdb_sha256"], "famdb_version_help": version_help}
        exported = outdir / "dfam.selected.hmm"
        manifest = outdir / "dfam.export.manifest.json"
        try:
            cached = json.loads(manifest.read_text()) if manifest.exists() else {}
        except (OSError, json.JSONDecodeError):
            cached = {}
        if (cached.get("signature") != signature or not exported.exists()
                or cached.get("sha256") != sha256(exported)):
            partial = Path(str(exported) + ".partial")
            with open(partial, "w") as handle:
                result = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, text=True, check=True)
            (outdir / "dfam.export.log").write_text(result.stderr)
            # FamDB may report missing partitions or unknown taxa yet exit 0.
            if result.stderr.strip():
                raise ValueError(f"FamDB export reported diagnostics; check dfam.export.log: {result.stderr.strip()}")
            if not parse_dfam_hmm_metadata(partial):
                raise FileNotFoundError("FamDB exported no nucleotide HMM profiles for the requested taxon")
            os.replace(partial, exported)
            write_json(manifest, {"signature": signature, "sha256": sha256(exported)})
        source = exported
        provenance.update(export_signature=signature, export_command=cmd)
    elif cfg.get("taxon"):
        raise ValueError("Taxonomic filtering requires a FamDB directory; use taxon: null for an already selected HMM export")
    if source.suffix == ".h5":
        raise ValueError("Pass the FamDB directory, not an individual HDF5 partition")
    metadata = parse_dfam_hmm_metadata(source)
    if not metadata:
        raise ValueError(f"No complete HMM profiles in {source}")
    if any(m["alphabet"] not in {"dna", "rna"} or m["model_length"] < 1 for m in metadata.values()):
        raise ValueError("Dfam requires nucleotide HMM profiles with ALPH DNA/RNA and LENG; protein HMMs are incompatible")
    provenance.update(hmm_sha256=sha256(source), n_models=len(metadata),
                      taxon=cfg["taxon"], curated_only=cfg["curated_only"])
    write_json(outdir / "dfam.reference.json", provenance)
    return source, metadata, provenance


def run_validation_search(config, fasta, outdir, logger, threads=1, search=True):
    cfg = options(config)
    if not cfg["enabled"]:
        return [], {"enabled": False, "policy_version": POLICY_VERSION}
    Path(outdir).mkdir(parents=True, exist_ok=True)
    try:
        reference, metadata, provenance = prepare_reference(cfg, outdir)
        nhmmer = executable(cfg.get("nhmmer_bin", "nhmmer"))
        version = capture([nhmmer, "-h"])
        if "nhmmer" not in version.lower():
            raise ValueError("Configured Dfam executable is not NHMMER")
        # A content digest plus executable bytes make the cache independent of
        # file timestamps and detect database replacement at the same path.
        db_size_mb = sum(len(record.seq) for record in SeqIO.parse(str(fasta), "fasta")) / 1e6
        if db_size_mb <= 0:
            raise ValueError("Dfam query FASTA contains no sequence")
        search_signature = {"db_size_mb": db_size_mb, "query_sha256": sha256(fasta), "hmm_sha256": provenance["hmm_sha256"],
                            "nhmmer_sha256": sha256(nhmmer), "nhmmer_help": version}
        search_dir = Path(outdir) / search_signature_digest(search_signature)
        tblout = search_dir / "dfam.merged.tblout"
        hits = []
        if search:
            tblout = run_nhmmer_parallel_chunks(reference, fasta, search_dir,
                        chunks=cfg.get("chunks", 1), nhmmer_bin=nhmmer,
                        cpu_per_job=cfg.get("threads", threads), max_parallel=cfg.get("jobs", 1),
                        db_size_mb=db_size_mb)
            # Keep weak matches in the evidence table; acceptance is decided later.
            hits = parse_nhmmer_tblout(tblout, max_evalue=float("inf"), min_score=-float("inf"), dfam_metadata=metadata)
            if (sha256(reference) != search_signature["hmm_sha256"]
                    or sha256(fasta) != search_signature["query_sha256"]
                    or sha256(nhmmer) != search_signature["nhmmer_sha256"]):
                raise RuntimeError("Dfam search inputs changed during execution; no completion was recorded")
        report = {"enabled": True, "status": "complete", "reference": provenance,
                  "search": search_signature, "tblout": str(tblout), "config": cfg,
                  "policy_version": POLICY_VERSION}
        if search:
            write_json(Path(outdir) / "dfam.validation.manifest.json", report)
        return hits, report
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        detail = str(exc)
        if isinstance(exc, subprocess.CalledProcessError):
            detail += ": " + (exc.stderr or exc.output or "no diagnostic output").strip()
            (Path(outdir) / "dfam.failure.log").write_text(detail + "\n")
        if cfg["on_missing"] == "error":
            raise
        logger.warning("Dfam validation unavailable: %s", detail)
        report = {"enabled": True, "status": "unavailable", "reason": detail,
                  "config": cfg, "policy_version": POLICY_VERSION}
        write_json(Path(outdir) / "dfam.validation.manifest.json", report)
        return [], report


def search_signature_digest(signature):
    return hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:20]


def rejection_reason(hit, cfg):
    if (not all(math.isfinite(v) for v in (hit.score, hit.evalue)) or hit.evalue < 0
            or (hit.threshold is not None and not math.isfinite(hit.threshold))):
        return "invalid_statistics"
    if cfg["curated_only"] and hit.curated is not True:
        return "not_curated"
    if hit.evalue > cfg["max_evalue"] or hit.score < cfg["min_score"]:
        return "weak_statistics"
    if hit.threshold is None and cfg["require_ga"]:
        return "missing_model_threshold"
    if hit.threshold is not None and hit.score < hit.threshold:
        return "below_model_threshold"
    if min(hit.hmm_start, hit.hmm_end, hit.ali_start, hit.ali_end) < 1:
        return "invalid_coordinates"
    if max(hit.hmm_start, hit.hmm_end) > hit.model_length or max(hit.ali_start, hit.ali_end) > hit.query_length:
        return "invalid_coordinates"
    if hit.query_coverage < cfg["min_query_coverage"] or hit.model_coverage < cfg["min_model_coverage"]:
        return "insufficient_coverage"
    if hit.query_coverage > 1 or hit.model_coverage > 1:
        return "invalid_coordinates"
    if hit.dfam_class not in {"Class_I", "Class_II"} or hit.dfam_order not in {"LTR", "LINE", "Penelope", "TIR", "Helitron", "SINE"}:
        return "unmapped_classification"
    return ""


def disagree(left, right):
    aliases = {"bel-pao": "pao", "gypsy-erv": "gypsy", "piggybac": "piggybac"}
    def normalized(key, value):
        if key != "superfamily":
            return value
        lower = str(value).lower()
        return aliases.get(lower, lower)
    return any(is_known(left[k]) and is_known(right[k])
               and normalized(k, left[k]) != normalized(k, right[k])
               for k in ("class", "order", "superfamily"))


def hit_label(hit):
    return dict(zip(("class", "order", "superfamily"),
                    (hit.dfam_class, hit.dfam_order, hit.dfam_superfamily)))


def integrate_validation(families, results, hits, config):
    """Apply reliable, unambiguous Dfam evidence after the primary classifier.

    Baselines must be built without Dfam. Existing known hierarchy nodes always
    survive contradictions. Confidence still comes from classify_family; Dfam
    alone is provisional, and correlated homology cannot promote confidence.
    """
    cfg = options(config)
    if not cfg["enabled"]:
        return results, []
    by_id = {f.family_id: f for f in families}
    grouped = {}
    evidence = []
    for hit in hits:
        reason = rejection_reason(hit, cfg)
        row = {**asdict(hit), "query_coverage": hit.query_coverage,
               "model_coverage": hit.model_coverage, "rejection_reason": reason,
               "accepted": not reason, "agreement": "not_used", "contributed": False}
        evidence.append(row)
        if not reason and hit.family_id in by_id:
            grouped.setdefault(hit.family_id, []).append((hit, row))
    final = []
    for baseline in results:
        row = dict(baseline)
        row.update({f"original_{k}": baseline[k] for k in ("class", "order", "superfamily", "confidence", "score", "evidence")})
        row.update(dfam_agreement="no_reliable_match", dfam_contributed=False)
        candidates = sorted(grouped.get(row["family_id"], []), key=lambda pair: (-pair[0].score, pair[0].evalue, pair[0].accession))
        if candidates:
            best = candidates[0][0]
            competing = [pair for pair in candidates if pair[0].score >= best.score * cfg["competition_score_fraction"]]
            ambiguous = any(disagree(hit_label(best), hit_label(pair[0])) for pair in competing)
            label = hit_label(best)
            conflict = disagree(baseline, label)
            agreement = "ambiguous" if ambiguous else "conflict" if conflict else "consistent"
            row["dfam_agreement"] = agreement
            if not ambiguous and not conflict:
                family = by_id[row["family_id"]]
                enriched = replace(family, dfam_class=best.dfam_class, dfam_order=best.dfam_order,
                                   dfam_superfamily=best.dfam_superfamily, dfam_model=best.model_name, dfam_score=best.score)
                proposed = classify_family(enriched)
                # Only fill missing hierarchy nodes; never erase a known node.
                for key in ("class", "order", "superfamily"):
                    if not is_known(row[key]) and is_known(label[key]):
                        row[key] = label[key]
                changed = any(row[k] != baseline[k] for k in ("class", "order", "superfamily"))
                if changed:
                    row["dfam_contributed"] = True
                    row["evidence"] += f";Dfam_validation={best.accession}:{best.model_name}"
                    if baseline["class"] == "Unknown":
                        row.update(status="putative", confidence="LOW" if proposed["confidence"] == "LOW" else "MEDIUM",
                                   score=min(0.64, proposed["score"]))
                # Explicit independence is needed, and the candidate must agree
                # with the final assignment, before using its confidence score.
                if (cfg["independent"] and not disagree(row, proposed)
                        and proposed["score"] > row["score"]):
                    row.update(score=proposed["score"], confidence=proposed["confidence"], status=proposed["status"])
                    row["dfam_contributed"] = True
                    if baseline["class"] == "Unknown":
                        row.update(status="putative", confidence="MEDIUM", score=min(0.64, row["score"]))
            if (row["dfam_contributed"]
                    and f"Dfam_validation={best.accession}:{best.model_name}" not in row["evidence"]):
                row["evidence"] += f";Dfam_validation={best.accession}:{best.model_name}"
            if agreement in {"ambiguous", "conflict"}:
                row["evidence"] += f";Dfam_validation_{agreement}"
            for hit, ev in candidates:
                ev["agreement"] = agreement if (hit, ev) in competing else "lower_scoring"
                ev["contributed"] = row["dfam_contributed"] and hit is best
        final.append(row)
    final_by_id = {r["family_id"]: r for r in final}
    for ev in evidence:
        result = final_by_id.get(ev["family_id"])
        if result:
            for key in ("class", "order", "superfamily", "confidence"):
                ev[f"original_{key}"] = result[f"original_{key}"]
                ev[f"final_{key}"] = result[key]
    return final, evidence


def write_validation_reports(evidence, outdir, prefix):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fields = list(DfamHit.__dataclass_fields__) + [
        "query_coverage", "model_coverage", "rejection_reason", "accepted", "agreement", "contributed",
        "original_class", "final_class", "original_order", "final_order",
        "original_superfamily", "final_superfamily", "original_confidence", "final_confidence",
    ]
    paths = {}
    for kind, rows in (("dfam_evidence", evidence), ("dfam_conflicts", [r for r in evidence if r["agreement"] in {"conflict", "ambiguous"}])):
        path = outdir / f"{prefix}.{kind}.tsv"
        with open(path, "w") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        paths[kind + "_tsv"] = str(path)
    return paths
