# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------

def val(f, name, default=0.0):
    return getattr(f, name, default)


def clamp(x):
    return max(0.0, min(float(x), 0.99))


def has_penelope_domain(domains):
    return any("penelope" in str(d).lower() for d in domains)


# ---------------------------------------------------------------------
# LTR
# ---------------------------------------------------------------------

def score_ltr(f):
    # --- STRUCTURE ---
    structure = 0.0

    if val(f, "ltr_present"):
        structure += 0.5

    ltr_type = val(f, "ltr_structural_type", "none")
    if ltr_type == "LTR_like":
        structure += 0.2
    elif ltr_type in {"LARD_like", "TRIM"}:
        structure += 0.1

    if val(f, "tg_ca_motif"):
        structure += 0.1

    if val(f, "ppt_like"):
        structure += 0.1

    structure = min(structure, 1.0)

    # --- DOMAIN ---
    domain = 0.0
    if val(f, "rt_present"):
        domain += 0.4
    if val(f, "integrase_present"):
        domain += 0.3
    if val(f, "rnaseh_present"):
        domain += 0.2
    if val(f, "gag_present"):
        domain += 0.1

    domain = min(domain, 1.0)

    # --- HOMOLOGY ---
    homology = 0.0
    if val(f, "homology_order") == "LTR":
        homology = val(f, "homology_score")
    elif val(f, "dfam_order") == "LTR":
        homology = min(1.0, val(f, "dfam_score") / 100.0)

    # --- FINAL ---
    score = (
        0.4 * structure +
        0.3 * domain +
        0.2 * homology +
        0.1 * val(f, "boundary_consistency")
    )

    # penalties
    if val(f, "tir_present") and not val(f, "ltr_present"):
        score -= 0.15

    if val(f, "dfam_order") in {"LINE", "Penelope"}:
        score -= 0.15

    # Canonical autonomous LTR support:
    # LTR structure + RT + integrase is strong evidence.
    if (
        val(f, "ltr_present")
        and val(f, "rt_present")
        and val(f, "integrase_present")
    ):
        score = max(score, 0.75)

    # Strong homology + LTR structure should also pass as OK.
    if (
        val(f, "ltr_present")
        and val(f, "homology_order") == "LTR"
        and val(f, "homology_score") >= 0.8
    ):
        score = max(score, 0.75)

    return clamp(score)


# ---------------------------------------------------------------------
# DNA TIR
# ---------------------------------------------------------------------

def _tsd_bonus(f):
    if not val(f, "tsd_present"):
        return 0.0

    tsd_seq = (val(f, "tsd_seq", "") or "").upper()
    tsd_len = int(val(f, "tsd_len", 0))
    support = float(val(f, "tsd_support", 0.0))

    bonus = 0.1

    if tsd_seq in {"TA", "TTAA", "TAA", "TTA"} or tsd_len in {2, 3, 8, 9, 10, 11}:
        bonus += 0.1

    if support >= 0.3:
        bonus += 0.05
    if support >= 0.5:
        bonus += 0.05

    return min(bonus, 0.25)


def score_dna_tir(f):
    # --- STRUCTURE ---
    structure = 0.0

    if val(f, "tir_present"):
        structure += 0.5

    if val(f, "tir_grade") in {"CREDIBLE", "STRONG"}:
        structure += 0.2

    structure += _tsd_bonus(f)
    structure = min(structure, 1.0)

    # --- DOMAIN ---
    domain = 0.0

    if val(f, "transposase_present"):
        domain += 0.8

    domains = val(f, "domains", set())
    specific = {"PIGGYBAC", "MUTATOR", "PIF_HARBINGER", "TCMAR", "HAT", "CACTA"}
    if domains & specific:
        domain += 0.2

    domain = min(domain, 1.0)

    # --- HOMOLOGY ---
    homology = 0.0
    if val(f, "homology_order") in {"DNA", "TIR"}:
        homology = val(f, "homology_score")
    elif val(f, "dfam_order") == "TIR":
        homology = min(1.0, val(f, "dfam_score") / 100.0)

    # --- FINAL ---
    score = (
        0.4 * structure +
        0.3 * domain +
        0.2 * homology +
        0.1 * val(f, "boundary_consistency")
    )

    # penalties
    if val(f, "rt_present") and not val(f, "transposase_present"):
        score -= 0.2

    return clamp(score)


# ---------------------------------------------------------------------
# LINE
# ---------------------------------------------------------------------

def score_line(f):
    # --- DOMAIN ---
    domain = 0.0

    if val(f, "rt_present"):
        domain += 0.7

    if val(f, "rnaseh_present"):
        domain += 0.2

    if val(f, "integrase_present"):
        domain -= 0.3

    if val(f, "ltr_present"):
        domain -= 0.2

    if has_penelope_domain(val(f, "domains", set())):
        domain -= 0.3

    domain = clamp(domain)

    # --- HOMOLOGY ---
    homology = 0.0
    if val(f, "homology_order") == "LINE":
        homology = val(f, "homology_score")
    elif val(f, "dfam_order") == "LINE":
        homology = min(1.0, val(f, "dfam_score") / 100.0)

    # --- FINAL ---
    score = (
        0.5 * domain +
        0.3 * homology +
        0.2 * (1.0 if val(f, "polyA_present") else 0.0)
    )

    if val(f, "tir_present"):
        score -= 0.2

    return clamp(score)


# ---------------------------------------------------------------------
# PENELOPE
# ---------------------------------------------------------------------

def score_penelope(f):
    domains = val(f, "domains", set())

    # --- DOMAIN ---
    if has_penelope_domain(domains):
        domain = 1.0
    elif val(f, "rt_present"):
        domain = 0.5
    else:
        domain = 0.0

    # --- HOMOLOGY ---
    homology = 0.0

    for v in [
        val(f, "homology_superfamily"),
        val(f, "dfam_superfamily"),
        val(f, "dfam_model"),
        val(f, "header_superfamily"),
    ]:
        if v and "penelope" in str(v).lower():
            homology = max(homology, 0.85)

    if val(f, "dfam_score") and homology > 0:
        homology = max(homology, min(0.99, val(f, "dfam_score") / 100.0))

    # --- PENALTIES ---
    penalty = 0.0

    if val(f, "ltr_present"):
        penalty += 0.25
    if val(f, "integrase_present"):
        penalty += 0.25
    if val(f, "tir_present"):
        penalty += 0.15

    # --- FINAL ---
    score = 0.4 * domain + 0.6 * homology - penalty

    return clamp(score)


# ---------------------------------------------------------------------
# HELITRON
# ---------------------------------------------------------------------

def score_helitron(f):
    structure = 1.0 if val(f, "helitron_signal") else 0.0

    domain = 0.8 if val(f, "helitron_domain_present") else 0.0

    homology = 0.0
    if val(f, "homology_order") == "Helitron":
        homology = val(f, "homology_score")
    elif val(f, "dfam_order") == "Helitron":
        homology = min(1.0, val(f, "dfam_score") / 100.0)

    score = (
        0.4 * structure +
        0.3 * domain +
        0.3 * homology
    )

    return clamp(score)


# ---------------------------------------------------------------------
# SINE
# ---------------------------------------------------------------------

def score_sine(f):
    homology = 0.0
    if val(f, "homology_order") == "SINE":
        homology = val(f, "homology_score")
    elif val(f, "dfam_order") == "SINE":
        homology = min(1.0, val(f, "dfam_score") / 100.0)

    polyA = 1.0 if val(f, "polyA_present") else 0.0
    poliii = min(1.0, val(f, "poliii_score"))
    sine_struct = min(1.0, val(f, "sine_score"))

    score = (
        0.5 * homology +
        0.2 * polyA +
        0.15 * poliii +
        0.15 * sine_struct
    )

    if val(f, "tir_present"):
        score -= 0.3

    if val(f, "ltr_present"):
        score -= 0.2

    return clamp(score)
