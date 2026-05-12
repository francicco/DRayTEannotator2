from dataclasses import dataclass, field
from typing import Set


RT_DOMAINS = {"RT"}

INTEGRASE_DOMAINS = {"INTEGRASE"}

TRANSPOSASE_DOMAINS = {"TRANSPOSASE"}

RH_DOMAINS = {"RNASEH"}

GAG_DOMAINS = {"GAG"}

HELITRON_DOMAINS = {"HELITRON_REP", "HELICASE"}


@dataclass
class Family:
    family_id: str
    consensus_len: int
    n_copies: int

    original_header: str = ""
    header_class: str = "Unknown"
    header_superfamily: str = "Unknown"

    homology_class: str = "Unknown"
    homology_order: str = "Unknown"
    homology_superfamily: str = "Unknown"
    homology_score: float = 0.0

    dfam_class: str | None = None
    dfam_order: str | None = None
    dfam_superfamily: str | None = None
    dfam_model: str | None = None
    dfam_score: float = 0.0

    rescue_class: str = "Unknown"
    rescue_order: str = "Unknown"
    rescue_superfamily: str = "Unknown"
    rescue_identity: float = 0.0
    rescue_aln_len: int = 0
    rescue_bits: float = 0.0
    rescue_target: str = ""

    domains: Set[str] = field(default_factory=set)

    ltr_present: bool = False
    tir_present: bool = False
    tir_len: int = 0
    tir_identity: float = 0.0
    tir_confidence: str = "LOW"
    tir_5p_motif: str = "NA"
    tir_3p_motif: str = "NA"
    tir_terminal_motif: str = "NA"
    tir_grade: str = "ABSENT"
    structural_type: str = "unknown"
    tir_tsd_element: bool = False
    structural_superfamily_hint: str = "unknown"
    structural_superfamily_confidence: str = "NONE"
    helitron_signal: bool = False

    tsd_present: bool = False
    tsd_seq: str = ""
    tsd_len: int = 0
    tsd_support: float = 0.0
    tsd_consensus: str = ""
    tsd_identity: float = 0.0
    polyA_present: bool = False
    polyA_tail_len: int = 0
    polyA_fraction: float = 0.0
    poliii_a_box: bool = False
    poliii_b_box: bool = False
    poliii_score: float = 0.0
    sine_score: float = 0.0
    sine_candidate: bool = False
    sine_confidence: str = "LOW"

    orf_count: int = 0
    orf_max_len: int = 0

    boundary_consistency: float = 0.0
    fragmentation_score: float = 1.0

    @property
    def rt_present(self) -> bool:
        return bool(self.domains & RT_DOMAINS)

    @property
    def integrase_present(self) -> bool:
        return bool(self.domains & INTEGRASE_DOMAINS)

    @property
    def transposase_present(self) -> bool:
        return bool(self.domains & TRANSPOSASE_DOMAINS)

    @property
    def rnaseh_present(self) -> bool:
        return bool(self.domains & RH_DOMAINS)

    @property
    def gag_present(self) -> bool:
        return bool(self.domains & GAG_DOMAINS)

    @property
    def helitron_domain_present(self) -> bool:
        return bool(self.domains & HELITRON_DOMAINS)
