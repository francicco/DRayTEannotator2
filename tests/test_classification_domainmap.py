from drayte.classification.domainmap import (
    normalize_domain_name,
    normalize_domains,
)


def test_normalize_domain_name():
    assert normalize_domain_name("RVT_1") == "RT"
    assert normalize_domain_name("Reverse_transcriptase") == "RT"
    assert normalize_domain_name("rve") == "INTEGRASE"
    assert normalize_domain_name("DDE_Tnp_1") == "TRANSPOSASE"
    assert normalize_domain_name("RNase_H") == "RNASEH"
    assert normalize_domain_name("unknown_domain") is None


def test_normalize_domains():
    domains = normalize_domains(
        ["RVT_1", "rve", "DDE_Tnp_1", "unknown_domain"]
    )

    assert domains == {"RT", "INTEGRASE", "TRANSPOSASE"}
