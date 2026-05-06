def is_ltr_candidate(f):
    return (
        f.ltr_present and
        (
            f.rt_present or
            f.integrase_present or
            f.homology_class == "LTR"
        )
    )


def is_dna_tir_candidate(f):
    return (
        f.tir_present or
        f.transposase_present
    )


def is_line_candidate(f):
    return (
        (
            f.homology_class == "LINE"
        )
        or
        (
            f.rt_present and
            not f.ltr_present and
            (
                f.polyA_present or
                f.homology_class == "LINE"
            )
        )
    )


def is_helitron_candidate(f):
    return (
        f.helitron_signal or
        f.homology_class == "Helitron"
    )
