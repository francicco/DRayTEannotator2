from dataclasses import dataclass

@dataclass
class Family:
    family_id: str
    consensus_len: int
    n_copies: int

    homology_class: str
    homology_superfamily: str
    homology_score: float

    rt_present: bool
    integrase_present: bool
    transposase_present: bool

    ltr_present: bool
    tir_present: bool
    helitron_signal: bool

    tsd_present: bool
    polyA_present: bool

    orf_count: int
    orf_max_len: int

    boundary_consistency: float
    fragmentation_score: float
