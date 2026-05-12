def is_ltr_candidate(f):
    if getattr(f, "ltr_present", False) and getattr(f, "ltr_structural_type", "none") in {
        "LTR_like",
        "LARD_like",
        "TRIM",
    }:
        return True

    if getattr(f, "rt_present", False):
        return True

    if getattr(f, "integrase_present", False):
        return True

    if getattr(f, "homology_class", "Unknown") == "LTR":
        return True

    if getattr(f, "homology_order", "Unknown") == "LTR":
        return True

    if getattr(f, "dfam_order", "Unknown") == "LTR":
        return True

    return False

def is_line_candidate(f):

    # RT but no LTR hallmarks
    if (
        f.rt_present and
        not f.ltr_present and
        not f.integrase_present
    ):
        return True

    # PolyA + RT-like
    if f.polyA_present and f.rt_present:
        return True

    # Homology
    if f.homology_class == "LINE":
        return True

    return False


def is_dna_tir_candidate(f):

    # Structural TIR
    if f.tir_present:
        return True

    # DDE transposase
    if f.transposase_present:
        return True

    # Homology
    if f.homology_class == "DNA":
        return True

    return False


def is_helitron_candidate(f):

    if f.helitron_signal:
        return True

    if f.homology_class == "Helitron":
        return True

    return False


def is_sine_candidate(f):

    if f.homology_class == "SINE":
        return True

    return False
