from drayte.pipeline.classification import (
    classification_input_signature,
    classification_signatures_compatible,
    classified_library_matches_input,
)


def test_classified_library_matches_input_ignores_classification_suffix(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">family-1#Unknown\nACGT\n")

    classified = tmp_path / "classified.fa"
    classified.write_text(">family-1#DNA/TcMar\nACGT\n")

    assert classified_library_matches_input(source, classified)


def test_classified_library_does_not_match_changed_sequence(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">family-1#Unknown\nACGT\n")

    classified = tmp_path / "classified.fa"
    classified.write_text(">family-1#DNA/TcMar\nAAAA\n")

    assert not classified_library_matches_input(source, classified)


def make_signature(tmp_path, classification_cfg):
    library = tmp_path / "library.fa"
    structure_library = tmp_path / "structure.fa"
    extensionwork = tmp_path / "extensionwork"
    extensionwork.mkdir(exist_ok=True)
    library.write_text(">family-1#Unknown\nACGT\n")
    structure_library.write_text(">family-1#Unknown\nACGT\n")
    return classification_input_signature(
        library=library,
        structure_library=structure_library,
        extensionwork_dir=extensionwork,
        classification_cfg=classification_cfg,
        stage_name="classification",
        final_mode=True,
    )


def test_execution_only_settings_do_not_invalidate_classification(tmp_path):
    cached = make_signature(
        tmp_path,
        {
            "min_orf_nt": 300,
            "chunks": 40,
            "jobs": 40,
            "cpu": 1,
        },
    )
    current = make_signature(
        tmp_path,
        {
            "min_orf_nt": 300,
            "chunks": 10,
            "jobs": 10,
            "cpu": 4,
        },
    )

    compatible, legacy_migration = classification_signatures_compatible(
        cached,
        current,
    )

    assert compatible
    assert not legacy_migration
    assert cached["execution_config"] != current["execution_config"]


def test_result_affecting_setting_invalidates_classification(tmp_path):
    cached = make_signature(tmp_path, {"min_orf_nt": 300, "chunks": 40})
    current = make_signature(tmp_path, {"min_orf_nt": 600, "chunks": 10})

    compatible, legacy_migration = classification_signatures_compatible(
        cached,
        current,
    )

    assert not compatible
    assert not legacy_migration


def test_legacy_signature_with_same_biological_inputs_can_be_migrated(tmp_path):
    current = make_signature(tmp_path, {"min_orf_nt": 300, "chunks": 10})
    legacy = {
        "schema_version": 1,
        **{key: current[key] for key in (
            "stage_name",
            "final_mode",
            "library_sha256",
            "structure_library_sha256",
            "extensionwork_sha256",
        )},
        "classification_config_sha256": "legacy-full-config-hash",
    }

    compatible, legacy_migration = classification_signatures_compatible(
        legacy,
        current,
    )

    assert compatible
    assert legacy_migration
