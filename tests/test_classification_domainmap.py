from drayte.classification.domainmap import (
    infer_orders_from_domains,
    infer_superfamilies_from_domains,
    normalize_domains,
)


def test_infer_orders_from_domains():
    domains = {"RVT_1", "rve", "DDE_Tnp_1", "RepHel"}

    orders = infer_orders_from_domains(domains)

    assert "RETRO" in orders
    assert "LTR" in orders
    assert "TIR" in orders
    assert "Helitron" in orders


def test_infer_superfamilies_from_domains():
    domains = {"PiggyBac", "MULE", "RepHel"}

    superfamilies = infer_superfamilies_from_domains(domains)

    assert "PiggyBac" in superfamilies
    assert "Mutator" in superfamilies
    assert "Helitron" in superfamilies
