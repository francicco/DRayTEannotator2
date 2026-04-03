def score_ltr(f):
    structure = 1.0 if f.ltr_present else 0.0

    domain = 0.0
    if f.rt_present:
        domain += 0.6
    if f.integrase_present:
        domain += 0.4

    homology = f.homology_score if f.homology_class == "LTR" else 0.0

    boundary = f.boundary_consistency
    annotation = 1.0 - f.fragmentation_score

    return (
        0.30 * structure +
        0.25 * domain +
        0.20 * homology +
        0.15 * boundary +
        0.10 * annotation
    )


def score_dna_tir(f):
    structure = 1.0 if f.tir_present else 0.0

    domain = 1.0 if f.transposase_present else 0.0

    homology = f.homology_score if f.homology_class == "DNA" else 0.0

    boundary = f.boundary_consistency
    annotation = 1.0 - f.fragmentation_score

    return (
        0.30 * structure +
        0.30 * domain +
        0.20 * homology +
        0.10 * boundary +
        0.10 * annotation
    )


def score_line(f):
    domain = 1.0 if f.rt_present else 0.0

    homology = f.homology_score if f.homology_class == "LINE" else 0.0

    polyA = 1.0 if f.polyA_present else 0.0

    boundary = f.boundary_consistency

    return (
        0.35 * domain +
        0.25 * homology +
        0.20 * polyA +
        0.20 * boundary
    )
