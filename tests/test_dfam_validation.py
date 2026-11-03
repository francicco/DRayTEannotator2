import csv
import json
import logging
from dataclasses import replace
from pathlib import Path

import pytest

from drayte.classification import dfam_validation as validation
from drayte.classification.dfam import DfamHit, parse_nhmmer_tblout
from drayte.classification.dfammap import parse_dfam_hmm_metadata
from drayte.classification.io import write_classification_tsv, write_evidence_tsv
from drayte.classification.models import Family
from drayte.classification.run import classify_families
from drayte.classification.write_library import rewrite_fasta_headers


def model(path, accession="DF000000001.1", alphabet="DNA", ga="30.0;", order="LINE", sf="CR1"):
    path.write_text(f"HMMER3/f\nNAME example\nACC {accession}\nLENG 100\nALPH {alphabet}\nGA {ga}\nCT Transposable_Element;Class_I_Retrotransposition;Non_LTR_Retrotransposon;CR1\nCC RepeatMasker Annotations:\nCC      Type: {order}\nCC      SubType: {sf}\n//\n")
    return path


def hit(**changes):
    base = DfamHit("fam", "example", "DF000000001.1", 1e-20, 60, 100, 1,
                   dfam_class="Class_I", dfam_order="LINE", dfam_superfamily="CR1",
                   hmm_start=1, hmm_end=100, query_length=200, model_length=100,
                   curated=True, threshold=30, threshold_source="GA")
    return replace(base, **changes)


def family(**changes):
    return Family("fam", 200, 0, **changes)


def integrate(f=None, hits=None, **cfg):
    families = [f or family()]
    baseline = classify_families(families)
    final, evidence = validation.integrate_validation(families, baseline, hits if hits is not None else [hit()], {"enabled": True, **cfg})
    return baseline[0], final[0], evidence


def test_real_nhmmer_column_orientation_and_metadata(tmp_path):
    metadata = parse_dfam_hmm_metadata(model(tmp_path / "reference.hmm"))
    tbl = tmp_path / "hits.tblout"
    tbl.write_text("fam#Unknown - example DF000000001.1 1 100 150 51 160 41 200 - 1e-20 60 0 description\n")
    parsed = parse_nhmmer_tblout(tbl, dfam_metadata=metadata)[0]
    assert parsed.family_id == "fam"
    assert parsed.query_coverage == .5
    assert parsed.model_coverage == 1
    assert parsed.threshold == 30
    assert parsed.curated is True
    assert parsed.dfam_order == "LINE"


@pytest.mark.parametrize("order,sf,expected", [("DNA", "hAT", ("Class_II", "TIR", "hAT")), ("RC", "Helitron", ("Class_II", "Helitron", "Helitron")), ("LINE", "Penelope", ("Class_I", "Penelope", "Penelope")), ("LTR", "Gypsy", ("Class_I", "LTR", "Gypsy")), ("SINE", "tRNA", ("Class_I", "SINE", "tRNA"))])
def test_authoritative_repeatmasker_metadata(tmp_path, order, sf, expected):
    meta = parse_dfam_hmm_metadata(model(tmp_path / "db.hmm", order=order, sf=sf))["example"]
    assert (meta["dfam_class"], meta["dfam_order"], meta["dfam_superfamily"]) == expected


def test_strong_unknown_assignment_is_provisional():
    baseline, final, evidence = integrate()
    assert baseline["class"] == final["original_class"] == "Unknown"
    assert (final["class"], final["order"], final["superfamily"]) == ("Class_I", "LINE", "CR1")
    assert final["confidence"] == "MEDIUM" and final["status"] == "putative"
    assert final["dfam_contributed"] and evidence[0]["contributed"]


@pytest.mark.parametrize("changes,reason", [({"evalue": 1}, "weak_statistics"), ({"score": 25}, "below_model_threshold"), ({"threshold": None}, "missing_model_threshold"), ({"hmm_end": 20}, "insufficient_coverage"), ({"query_length": 2000}, "insufficient_coverage"), ({"curated": False}, "not_curated"), ({"curated": None}, "not_curated"), ({"dfam_order": "Unknown"}, "unmapped_classification"), ({"score": float("nan")}, "invalid_statistics")])
def test_weak_hits_never_rescue(changes, reason):
    baseline, final, evidence = integrate(hits=[hit(**changes)])
    assert final["class"] == baseline["class"] == "Unknown"
    assert not final["dfam_contributed"]
    assert evidence[0]["rejection_reason"] == reason


def test_uncurated_can_be_enabled():
    _, final, _ = integrate(hits=[hit(curated=False)], curated_only=False)
    assert final["order"] == "LINE"


def test_no_reliable_hits_preserves_unknown():
    _, final, evidence = integrate(hits=[])
    assert final["class"] == "Unknown" and not evidence


def test_conflict_preserves_structural_domain_assignment():
    f = family(ltr_present=True, domains={"RT", "INTEGRASE"})
    baseline, final, evidence = integrate(f)
    for key in ("class", "order", "superfamily", "confidence", "score"):
        assert final[key] == baseline[key]
    assert final["dfam_agreement"] == "conflict"
    assert not final["dfam_contributed"] and evidence[0]["agreement"] == "conflict"


def test_competing_superfamilies_do_not_rescue():
    _, final, evidence = integrate(hits=[hit(), hit(accession="DF2", model_name="other", score=59, dfam_superfamily="RTE")])
    assert final["class"] == "Unknown" and final["dfam_agreement"] == "ambiguous"
    assert all(not ev["contributed"] for ev in evidence)


def test_noncompeting_hit_does_not_block_best():
    _, final, evidence = integrate(hits=[hit(), hit(accession="DF2", model_name="other", score=40, dfam_superfamily="RTE")])
    assert final["superfamily"] == "CR1" and evidence[1]["agreement"] == "lower_scoring"


def test_concordance_does_not_double_count_header_homology():
    f = family(homology_order="LINE", homology_superfamily="CR1", homology_score=.9)
    baseline, final, _ = integrate(f)
    assert final["confidence"] == baseline["confidence"]
    assert final["score"] == baseline["score"]
    assert not final["dfam_contributed"]


def test_order_and_class_level_refinement():
    families = [family()]
    for baseline in (
        {"family_id": "fam", "class": "Class_I", "order": "Unknown", "superfamily": "Unknown", "confidence": "LOW", "score": .4, "status": "putative", "evidence": "baseline"},
        {"family_id": "fam", "class": "Class_I", "order": "LINE", "superfamily": "Unknown", "confidence": "MEDIUM", "score": .6, "status": "likely", "evidence": "baseline"},
    ):
        final, _ = validation.integrate_validation(families, [baseline], [hit()], {"enabled": True})
        assert final[0]["order"] == "LINE" and final[0]["superfamily"] == "CR1"
        assert final[0]["original_evidence"] == "baseline"


def test_disabled_is_exact_baseline_and_skips_dependencies(tmp_path, monkeypatch):
    def unexpected(*args, **kwargs):
        pytest.fail("disabled validation must not access databases or tools")
    monkeypatch.setattr(validation, "prepare_reference", unexpected)
    assert validation.run_validation_search({"enabled": False}, "missing.fa", tmp_path, logging.getLogger())[0] == []
    families = [family()]
    baseline = classify_families(families)
    result, evidence = validation.integrate_validation(families, baseline, [hit()], {"enabled": False})
    assert result is baseline and evidence == []


def test_missing_reference_error_and_explicit_warning(tmp_path):
    with pytest.raises(FileNotFoundError, match="does not exist"):
        validation.run_validation_search({"enabled": True, "database": str(tmp_path / "missing")}, "missing.fa", tmp_path, logging.getLogger())
    hits, report = validation.run_validation_search({"enabled": True, "database": str(tmp_path / "missing"), "on_missing": "warn"}, "missing.fa", tmp_path, logging.getLogger())
    assert not hits and report["status"] == "unavailable"


def test_protein_database_rejected(tmp_path):
    with pytest.raises(ValueError, match="nucleotide HMM"):
        validation.prepare_reference(validation.options({"database": str(model(tmp_path / "protein.hmm", alphabet="amino"))}), tmp_path / "out")


def test_taxon_on_export_requires_explicit_preselection(tmp_path):
    with pytest.raises(ValueError, match="already selected"):
        validation.prepare_reference(validation.options({"database": str(model(tmp_path / "dna.hmm")), "taxon": "Arthropoda"}), tmp_path / "out")


def test_reports_and_repeatmasker_fasta_preserve_ids_and_sequences(tmp_path):
    f = family()
    _, final, evidence = integrate(f)
    table = tmp_path / "classifications.tsv"
    write_classification_tsv([final], table)
    write_evidence_tsv([f], [final], tmp_path / "evidence.tsv")
    paths = validation.write_validation_reports(evidence, tmp_path, "test")
    assert Path(paths["dfam_evidence_tsv"]).read_text().count("\n") == 2
    assert Path(paths["dfam_conflicts_tsv"]).read_text().count("\n") == 1
    fasta = tmp_path / "families.fa"
    fasta.write_text(">fam#Unknown\nACGTACGT\n")
    output = tmp_path / "classified.fa"
    rewrite_fasta_headers(fasta, table, output, keep_unknown=True)
    assert output.read_text().splitlines() == [">fam#LINE/CR1", "ACGTACGT"]


def test_nhmmer_search_reuse_and_database_content_invalidation(tmp_path, monkeypatch):
    db = model(tmp_path / "reference.hmm")
    fasta = tmp_path / "input.fa"
    fasta.write_text(">fam\n" + "A" * 200 + "\n")
    calls = []
    def fake_search(reference, source, outdir, **kwargs):
        calls.append(outdir)
        Path(outdir).mkdir(parents=True, exist_ok=True)
        tbl = Path(outdir) / "hits.tblout"
        tbl.write_text("fam - example DF000000001.1 1 100 1 100 1 100 200 + 1e-20 60 0 desc\n")
        return tbl
    monkeypatch.setattr(validation, "run_nhmmer_parallel_chunks", fake_search)
    monkeypatch.setattr(validation, "executable", lambda name: str(db))
    monkeypatch.setattr(validation, "capture", lambda cmd: "nhmmer 3.4")
    cfg = {"enabled": True, "database": str(db)}
    _, first = validation.run_validation_search(cfg, fasta, tmp_path / "out", logging.getLogger())
    validation.run_validation_search(cfg, fasta, tmp_path / "out", logging.getLogger())
    assert calls[0] == calls[1]
    db.write_text(db.read_text().replace("GA 30", "GA 35"))
    _, second = validation.run_validation_search(cfg, fasta, tmp_path / "out", logging.getLogger())
    assert calls[-1] != calls[0]
    assert first["reference"]["hmm_sha256"] != second["reference"]["hmm_sha256"]


def test_warning_mode_creates_new_output_directory(tmp_path):
    hits, report = validation.run_validation_search({"enabled": True, "on_missing": "warn"}, "missing.fa", tmp_path / "new", logging.getLogger())
    assert report["status"] == "unavailable" and not hits
    assert (tmp_path / "new" / "dfam.validation.manifest.json").exists()


@pytest.mark.parametrize("config", [{"min_query_coverage": 0}, {"min_model_coverage": float("nan")}, {"competition_score_fraction": 1.1}, {"engine": "hmmscan"}, {"enabled": "false"}, {"threads": 0}])
def test_invalid_options_are_actionable(config):
    with pytest.raises(ValueError, match="classification.dfam"):
        validation.options(config)


def test_pipeline_refreshes_only_classification_for_policy_and_database_changes(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from drayte.pipeline import classification as stage
    fasta = tmp_path / "library.fa"
    fasta.write_text(">fam#Unknown\n" + "ACGT" * 50 + "\n")
    cfg = SimpleNamespace(species="Test", outdir_path=tmp_path, threads=1,
                          extra={"classification": {"dfam": {"enabled": True}}})
    searches = []
    database = {"sha256": "first"}
    def fake_search(config, fasta, outdir, logger, threads=1, search=True):
        searches.append(dict(config))
        return [hit()], {"enabled": True, "status": "complete", "config": config.copy(), "reference": database.copy()}
    monkeypatch.setattr(stage, "run_validation_search", fake_search)
    monkeypatch.setattr(stage, "detect_tirs_from_fasta", lambda *a, **k: [])
    monkeypatch.setattr(stage, "detect_ltrs_from_fasta", lambda *a, **k: [])
    monkeypatch.setattr(stage, "detect_sines_from_fasta", lambda *a, **k: [])
    first = stage.run(cfg, {"final_library": str(fasta)}, logging.getLogger(), final_mode=True)
    assert first["n_classified"] == 1
    saved = Path(first["classifications_tsv"]).stat().st_mtime_ns
    def unexpected(*a, **k):
        pytest.fail("compatible classification checkpoint must skip rebuilding structures")
    monkeypatch.setattr(stage, "detect_complexity_from_fasta", unexpected)
    second = stage.run(cfg, {"final_library": str(fasta)}, logging.getLogger(), final_mode=True)
    assert Path(second["classifications_tsv"]).stat().st_mtime_ns == saved
    # Policy changes regenerate decisions while retaining structure evidence.
    monkeypatch.setattr(stage, "detect_complexity_from_fasta", lambda *a, **k: [])
    cfg.extra["classification"]["dfam"]["min_score"] = 100
    third = stage.run(cfg, {"final_library": str(fasta)}, logging.getLogger(), final_mode=True)
    assert third["n_classified"] == 0
    assert first["input_signature"]["dfam_validation_sha256"] != third["input_signature"]["dfam_validation_sha256"]
    database["sha256"] = "second"
    fourth = stage.run(cfg, {"final_library": str(fasta)}, logging.getLogger(), final_mode=True)
    assert fourth["input_signature"]["dfam_validation_sha256"] != third["input_signature"]["dfam_validation_sha256"]
    cfg.extra["classification"]["dfam"] = {"enabled": False}
    last = stage.run(cfg, {"final_library": str(fasta)}, logging.getLogger(), final_mode=True)
    assert last["n_classified"] == 0
    assert not last["dfam_validation"].get("status") == "unavailable"


def test_equivalent_superfamily_aliases_are_not_conflicts():
    f = family(homology_order="LTR", homology_superfamily="Bel-Pao", homology_score=.9)
    baseline, final, _ = integrate(f, hits=[hit(dfam_order="LTR", dfam_superfamily="Pao")])
    assert final["dfam_agreement"] == "consistent"
    assert final["superfamily"] == baseline["superfamily"] == "Bel-Pao"


def test_explicit_independence_reuses_classifier_confidence_and_records_contribution():
    f = family(homology_order="LINE", homology_superfamily="CR1", homology_score=.8)
    baseline, final, evidence = integrate(f, hits=[hit(score=95)], independent=True)
    assert final["score"] >= baseline["score"]
    assert final["dfam_contributed"] and evidence[0]["contributed"]
    assert "Dfam_validation=" in final["evidence"]
    assert final["original_evidence"] == baseline["evidence"]


def test_missing_nhmmer_has_actionable_status(tmp_path, monkeypatch):
    db = model(tmp_path / "db.hmm")
    def unavailable(name):
        raise FileNotFoundError("Dfam dependency unavailable: nhmmer")
    monkeypatch.setattr(validation, "executable", unavailable)
    _, report = validation.run_validation_search({"enabled": True, "database": str(db), "on_missing": "warn"}, "missing.fa", tmp_path / "out", logging.getLogger())
    assert report["status"] == "unavailable" and "nhmmer" in report["reason"]


def test_empty_profile_file_is_rejected(tmp_path):
    db = tmp_path / "empty.hmm"
    db.write_text("")
    with pytest.raises(ValueError, match="No complete HMM profiles"):
        validation.prepare_reference(validation.options({"database": str(db)}), tmp_path / "out")


def test_disabled_standalone_ignores_legacy_dfam_cache(tmp_path):
    import subprocess
    import sys
    fasta = tmp_path / "library.fa"
    fasta.write_text(">fam#Unknown\nACGTACGT\n")
    cfg = tmp_path / "dfam.yaml"
    cfg.write_text("enabled: false\n")
    stale = tmp_path / "stale"
    stale.mkdir()
    (stale / "dfam.merged.tblout").write_text("bad legacy evidence must not be parsed\n")
    output = tmp_path / "classifications.tsv"
    result = subprocess.run([sys.executable, "-m", "drayte.classification.run", "--fasta", str(fasta), "--dfam-config", str(cfg), "--dfam-outdir", str(stale), "--pfam-outdir", str(tmp_path / "pfam"), "--output", str(output)], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert rows[0]["class"] == "Unknown"
    assert "original_class" not in rows[0]


@pytest.mark.parametrize("order,cls,sf", [("SINE", "Class_I", "tRNA"), ("LTR", "Class_I", "Gypsy"), ("TIR", "Class_II", "hAT"), ("Helitron", "Class_II", "Helitron"), ("Penelope", "Class_I", "Penelope")])
def test_supported_orders_use_existing_scoring(order, cls, sf):
    _, final, _ = integrate(hits=[hit(dfam_order=order, dfam_class=cls, dfam_superfamily=sf)])
    assert final["order"] == order and final["class"] == cls
    assert .55 <= final["score"] <= .64
    assert final["status"] == "putative" and final["confidence"] == "MEDIUM"
