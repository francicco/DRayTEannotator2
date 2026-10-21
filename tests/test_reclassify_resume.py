from drayte.pipeline.reclassify import classified_output_matches_input


def test_classified_output_can_change_class_without_changing_family_sequence(tmp_path):
    merged = tmp_path / "merged.fa"
    merged.write_text(
        ">Pmex-rnd-1_family-0#Unknown\nACGT\n"
        ">Pmex-rnd-1_family-1#Unknown\nTTAA\n"
    )

    classified = tmp_path / "merged.fa.classified"
    classified.write_text(
        ">Pmex-rnd-1_family-0#DNA/TcMar\nACGT\n"
        ">Pmex-rnd-1_family-1#LINE/L1\nTTAA\n"
    )

    assert classified_output_matches_input(merged, classified)


def test_classified_output_is_rejected_when_a_sequence_changed(tmp_path):
    merged = tmp_path / "merged.fa"
    merged.write_text(">Pmex-rnd-1_family-0#Unknown\nACGT\n")

    classified = tmp_path / "merged.fa.classified"
    classified.write_text(">Pmex-rnd-1_family-0#DNA/TcMar\nAAAA\n")

    assert not classified_output_matches_input(merged, classified)
