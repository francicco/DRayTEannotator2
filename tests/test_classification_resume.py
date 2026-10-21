from drayte.pipeline.classification import classified_library_matches_input


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
