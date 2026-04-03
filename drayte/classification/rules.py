def is_ltr_candidate(f):
    return (
        f.ltr_present and
        (f.rt_present or f.homology_class == "LTR")
    )


def is_dna_tir_candidate(f):
    return (
        f.tir_present or
        f.transposase_present
    )


def is_line_candidate(f):
    return (
        f.rt_present or
        (f.polyA_present and f.homology_class == "LINE")
    )
