def score_ltr(f):
    structure = 0.0

    if f.ltr_present:
        structure += 0.45

    if getattr(f, "ltr_structural_type", "none") == "LTR_like":
        structure += 0.15
    elif getattr(f, "ltr_structural_type", "none") in {"LARD_like", "TRIM"}:
        structure += 0.10

    if getattr(f, "tg_ca_motif", False):
        structure += 0.10

    if getattr(f, "ppt_like", False):
        structure += 0.05

    structure = min(structure, 1.0)

    domain = 0.0
    if f.rt_present:
        domain += 0.45
    if f.integrase_present:
        domain += 0.25
    if getattr(f, "rnaseh_present", False):
        domain += 0.20
    if getattr(f, "gag_present", False):
        domain += 0.10

    domain = min(domain, 1.0)

    homology = f.homology_score if f.homology_class == "LTR" else 0.0

    boundary = f.boundary_consistency
    annotation = 1.0 - f.fragmentation_score

    return round(
        0.35 * structure +
        0.30 * domain +
        0.20 * homology +
        0.10 * boundary +
        0.05 * annotation,
        3,
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


def score_helitron(f):

    structure = 1.0 if f.helitron_signal else 0.0

    homology = (
        f.homology_score
        if f.homology_class == "Helitron"
        else 0.0
    )

    boundary = f.boundary_consistency

    annotation = 1.0 - f.fragmentation_score

    return (
        0.45 * structure +
        0.30 * homology +
        0.15 * boundary +
        0.10 * annotation
    )

def score_sine(f):
    homology = (
        f.homology_score
        if f.homology_class == "SINE"
        else 0.0
    )

    polyA = 1.0 if f.polyA_present else 0.0

    return (
        0.80 * homology +
        0.20 * polyA
    )

