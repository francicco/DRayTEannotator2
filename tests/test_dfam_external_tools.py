"""Opt-in export/search test against installed FamDB and NHMMER.

See docs/dfam-validation.md. No reference databases are downloaded by tests.
"""
import logging
import os
from pathlib import Path

import pytest

from drayte.classification.dfam_validation import (
    integrate_validation, run_validation_search,
)
from drayte.classification.features import build_families_from_evidence
from drayte.classification.run import classify_families


def test_real_famdb_export_nhmmer_and_resume(tmp_path):
    required = ("DRAYTE_FAMDB_TEST_DATABASE", "DRAYTE_FAMDB_BIN", "DRAYTE_NHMMER_BIN")
    if not all(os.environ.get(key) for key in required):
        pytest.skip("Set the documented DRAYTE_FAMDB_TEST_DATABASE/tool paths to run real FamDB/NHMMER integration")
    fasta = Path(__file__).parent / "data" / "dfam" / "consensuses.fa"
    cfg = {"enabled": True, "database": os.environ[required[0]],
           "famdb_bin": os.environ[required[1]], "nhmmer_bin": os.environ[required[2]],
           "taxon": os.environ.get("DRAYTE_FAMDB_TEST_TAXON", "1")}
    hits, report = run_validation_search(cfg, fasta, tmp_path, logging.getLogger())
    assert report["status"] == "complete" and hits
    families = build_families_from_evidence(fasta, domain_hits=[], dfam_hits=[])
    results, _ = integrate_validation(families, classify_families(families), hits, cfg)
    by_id = {row["family_id"]: row for row in results}
    assert by_id["MADE1_test"]["order"] == "TIR"
    assert by_id["MADE1_test"]["original_class"] == "Unknown"
    assert by_id["MADE1_test"]["status"] == "putative"
    assert by_id["unresolved_short"]["class"] == "Unknown"
    tblouts = list(tmp_path.glob("*/tblouts/*.tblout"))
    timestamps = {p: p.stat().st_mtime_ns for p in tblouts}
    rerun, second = run_validation_search(cfg, fasta, tmp_path, logging.getLogger())
    assert rerun == hits and report == second
    assert timestamps == {p: p.stat().st_mtime_ns for p in tblouts}
