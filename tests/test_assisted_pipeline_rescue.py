"""Exercise reference rescue through real feature/decision/header-writing code."""
import csv
import json
import logging
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio import SeqIO

from drayte.classification.structure_merge import StructureEvidence, write_structure_summary
from drayte.pipeline import classification as stage
from drayte.pipeline.config import PipelineConfig


@pytest.fixture
def assisted_case(tmp_path, monkeypatch):
    root = tmp_path / "T_DRayTEannotator2"
    outdir = root / "classification"
    outdir.mkdir(parents=True)
    sequence = "ACGTCATGCTAGGCTA" * 63
    other_sequence = "TGCATGACTAGCTCGA" * 63
    curated = root / "curated.fa"
    curated.write_text(f">T-family-1#Unknown\n{sequence}\n")
    structure = root / "broad_structure.fa"
    structure.write_text(
        f">T-family-1#Unknown\n{sequence}\n>T-family-2#Unknown\n{other_sequence}\n"
    )
    reference = tmp_path / "Heliconius.classification.fa"
    reference.write_text(f">Curated-Helitron#RC/Helitron\n{sequence}\n")
    calls_file = tmp_path / "mmseqs_calls.jsonl"
    fake_mmseqs = tmp_path / "fake_mmseqs"
    fake_mmseqs.write_text(
        f"#!{sys.executable}\n"
        "import json, sys\n"
        "from pathlib import Path\n"
        "def headers(path):\n"
        "    return [line[1:].split()[0] for line in Path(path).read_text().splitlines() if line.startswith('>')]\n"
        "query, target = headers(sys.argv[2]), headers(sys.argv[3])\n"
        "query_id = next((item for item in query if item.startswith('T-family-1#')), query[0])\n"
        "target_id = next((item for item in target if item.startswith('T-family-1#')), target[0])\n"
        "Path(sys.argv[4]).write_text(f'{query_id}\\t{target_id}\\t0.95\\t900\\t0\\t0\\t1\\t900\\t1\\t900\\t1e-80\\t1800\\n')\n"
        f"with Path({str(calls_file)!r}).open('a') as handle:\n"
        "    handle.write(json.dumps({'query': query, 'target': target}) + '\\n')\n"
    )
    fake_mmseqs.chmod(0o755)
    config = PipelineConfig(
        genome=str(tmp_path / "genome.fa"), outdir=str(root), species="T", threads=2,
        extra={"classification": {
            "structure_input_fasta": str(structure), "min_orf_nt": 300,
            "mmseqs_rescue_db": str(reference), "mmseqs_bin": str(fake_mmseqs),
        }},
    )
    extensionwork = root / "extension/extensionwork"
    extensionwork.mkdir(parents=True)
    structure_cache = outdir / "T.structure_summary.tsv"
    write_structure_summary([
        StructureEvidence("T-family-1", consensus_len=len(sequence)),
        StructureEvidence("T-family-2", consensus_len=len(other_sequence)),
    ], structure_cache)
    cached_paths = [structure_cache]
    for relative in ("pfam/domains.domtblout", "dfam/dfam.merged.tblout"):
        path = outdir / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# completed cached evidence: no hits\n")
        cached_paths.append(path)

    def unexpected(*args, **kwargs):
        pytest.fail("Existing completed structural/HMM evidence must be reused")

    monkeypatch.setattr(stage, "detect_complexity_from_fasta", lambda *args: [])
    for name in ("detect_tirs_from_fasta", "detect_sines_from_fasta", "detect_ltrs_from_fasta",
                 "run_domain_annotation", "run_nhmmer", "run_nhmmer_parallel_chunks"):
        monkeypatch.setattr(stage, name, unexpected)

    def signature():
        return stage.classification_input_signature(
            library=curated, structure_library=structure, extensionwork_dir=extensionwork,
            classification_cfg=config.extra["classification"], stage_name="classification",
            final_mode=True,
        )

    def run():
        return stage.run(config, {"final_library": str(curated)},
                         logging.getLogger("test.assisted.pipeline"), final_mode=True)

    def calls():
        return [json.loads(line) for line in calls_file.read_text().splitlines()] if calls_file.exists() else []

    return SimpleNamespace(root=root, outdir=outdir, curated=curated, structure=structure,
                           reference=reference, config=config, run=run, signature=signature,
                           cached_paths=cached_paths, calls=calls)


def legacy_outputs(case):
    """Reproduce the old completed-stage shape without rescue provenance."""
    signature = case.signature()
    signature.pop("mmseqs_reference_sha256", None)
    signature.pop("mmseqs_rescue_algorithm_version", None)
    signature.pop("mmseqs_rescue_mode", None)
    signature.pop("mmseqs_rescue_score_fraction", None)
    result = {"stage": "classification", "input_signature": signature,
              "input_library": str(case.curated), "n_classified": 0, "n_unknown": 2}
    for key, filename in (("structure_evidence", "T.structure_summary.tsv"),
                          ("classifications_tsv", "T.classifications.tsv"),
                          ("evidence_tsv", "T.evidence.tsv")):
        path = case.outdir / filename
        if not path.exists():
            path.write_text("legacy completed output\n")
        result[key] = str(path)
    classified = case.outdir / "T.classified.fa"
    classified.write_text(case.structure.read_text())
    result["classified_library"] = str(classified)
    (case.outdir / "classification.manifest.json").write_text(json.dumps(result))
    old_forward = case.outdir / "mmseqs/mmseqs_rescue.tsv"
    old_forward.parent.mkdir()
    old_forward.write_text(
        "T-family-1#Unknown\tCurated-Helitron#RC/Helitron\t0.95\t900\t0\t0\t1\t900\t1\t900\t1e-80\t1800\n"
    )
    return result


def read_table(path):
    with Path(path).open() as handle:
        return {row["family_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def test_old_assisted_result_is_refreshed_and_reference_labels_reach_final_fasta(assisted_case):
    case = assisted_case
    legacy_outputs(case)
    original_sequences = {record.id.split("#")[0]: str(record.seq)
                          for record in SeqIO.parse(case.structure, "fasta")}
    caches = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in case.cached_paths}
    result = case.run()

    assert len(case.calls()) == 1  # An unverified old TSV is never adopted.
    assert case.calls()[0]["query"] == ["T-family-1#Unknown", "T-family-2#Unknown"]
    assert case.calls()[0]["target"] == ["Curated-Helitron#RC/Helitron"]
    assert {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in case.cached_paths} == caches

    classifications = read_table(result["classifications_tsv"])
    assert set(classifications) == set(original_sequences)
    assert classifications["T-family-1"]["order"] == "Helitron"
    assert classifications["T-family-2"]["order"] == "Unknown"
    evidence = read_table(result["evidence_tsv"])
    assert evidence["T-family-1"]["rescue_order"] == "Helitron"
    assert evidence["T-family-1"]["rescue_target"] == "Curated-Helitron#RC/Helitron"
    output = list(SeqIO.parse(result["classified_library"], "fasta"))
    assert {record.id.split("#")[0]: str(record.seq) for record in output} == original_sequences
    assert next(record.id for record in output if record.id.startswith("T-family-1#")) == "T-family-1#RC/Helitron"
    assert result["n_classified"] == 1 and result["n_unknown"] == 1
    report = result["mmseqs_reference_rescue"]
    assert report["forward_hits"] == 1 and report["families_with_rescue"] == 1
    assert "reverse_hits" not in report and "reciprocal_pairs" not in report
    assert "mmseqs_reciprocal_rescue" not in result
    assert report["query_sha256"] == stage.sha256_file(case.structure)
    assert report["reference_sha256"] == stage.sha256_file(case.reference)
    assert result["mmseqs_rescue_impact"]["classifications_changed_by_rescue"] == 1
    assert result["mmseqs_rescue_impact"]["previously_unknown_rescued"] == 1


def test_matching_assisted_checkpoint_skips_whole_stage_after_worker_change(assisted_case, monkeypatch):
    case = assisted_case
    first = case.run()
    case.config.threads = 8
    case.config.extra["classification"].update(cpu=4, jobs=3, chunks=12)

    def unexpected(*args, **kwargs):
        pytest.fail("Valid assisted stage must skip all evidence work")

    monkeypatch.setattr(stage, "run_reference_mmseqs_rescue", unexpected)
    monkeypatch.setattr(stage, "detect_complexity_from_fasta", unexpected)
    second = case.run()
    assert second["classified_library"] == first["classified_library"]
    assert len(case.calls()) == 1


def test_reference_changed_in_place_invalidates_assisted_stage(assisted_case):
    case = assisted_case
    case.run()
    case.reference.write_text(case.reference.read_text().replace("ACGTC", "TCGTC", 1))
    with pytest.raises(RuntimeError, match="classification.mmseqs_rescue_db contents"):
        case.run()
    assert len(case.calls()) == 1


def test_changed_reference_regenerates_search_when_stage_checkpoint_is_absent(assisted_case):
    case = assisted_case
    case.run()
    caches = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in case.cached_paths}
    case.reference.write_text(case.reference.read_text().replace("ACGTC", "TCGTC", 1))
    (case.outdir / "classification.manifest.json").unlink()
    result = case.run()
    assert len(case.calls()) == 2
    assert {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in case.cached_paths} == caches
    report = result["mmseqs_reference_rescue"]
    assert report["reference_sha256"] == stage.sha256_file(case.reference)
    assert not report["cache_reused"]


def test_unassisted_checkpoint_still_skips_without_external_report(assisted_case, monkeypatch):
    case = assisted_case
    case.config.extra["classification"].pop("mmseqs_rescue_db")
    legacy_outputs(case)
    (case.outdir / "mmseqs/mmseqs_rescue.tsv").unlink()

    def unexpected(*args, **kwargs):
        pytest.fail("Unassisted completed stage should retain its old resume behavior")

    monkeypatch.setattr(stage, "run_reference_mmseqs_rescue", unexpected)
    monkeypatch.setattr(stage, "detect_complexity_from_fasta", unexpected)
    result = case.run()
    assert result["n_unknown"] == 2
    assert case.calls() == []
