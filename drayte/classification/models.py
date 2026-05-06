from dataclasses import dataclass, field
from typing import Set


RT_DOMAINS = {
    "RT",
    "RT_domain",
    "Reverse_transcriptase",
    "RVT_1",
    "RVT_2",
}

INTEGRASE_DOMAINS = {
    "INT",
    "INT_domain",
    "Integrase",
    "rve",
}

TRANSPOSASE_DOMAINS = {
    "Transposase",
    "Transposase_domain",
    "DDE_Tnp",
    "hAT",
    "Tc1",
    "Mutator",
}

RH_DOMAINS = {
    "RH",
    "RNaseH",
    "RNase_H",
}

GAG_DOMAINS = {
    "GAG",
    "Gag",
}

HELITRON_DOMAINS = {
    "Helitron",
    "RepHel",
    "Rep",
    "Helicase",
}


@dataclass
class Family:
    family_id: str
    consensus_len: int
    n_copies: int

    homology_class: str = "Unknown"
    homology_superfamily: str = "Unknown"
    homology_score: float = 0.0

    domains: Set[str] = field(default_factory=set)

    ltr_present: bool = False
    tir_present: bool = False
    helitron_signal: bool = False

    tsd_present: bool = False
    polyA_present: bool = False

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
