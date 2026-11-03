from __future__ import annotations

from pathlib import Path

def infer_te_from_dfam_model(
    model_name: str,
) -> tuple[str, str, str]:

    name = model_name.lower()

    if "helitron" in name:
        return "Class_II", "Helitron", "Helitron"

    if any(x in name for x in ["gypsy", "copia", "bel", "pao", "ltr"]):
        if "gypsy" in name:
            return "Class_I", "LTR", "Gypsy"
        if "copia" in name:
            return "Class_I", "LTR", "Copia"
        if "bel" in name or "pao" in name:
            return "Class_I", "LTR", "Pao"
        return "Class_I", "LTR", "Unknown"

    if any(x in name for x in ["line", "rte", "cr1", "l1", "r1", "r2", "jockey"]):
        if "cr1" in name:
            return "Class_I", "LINE", "CR1"
        if "rte" in name:
            return "Class_I", "LINE", "RTE"
        if "l1" in name:
            return "Class_I", "LINE", "L1"
        if "r1" in name:
            return "Class_I", "LINE", "R1"
        if "r2" in name:
            return "Class_I", "LINE", "R2"
        if "jockey" in name:
            return "Class_I", "LINE", "Jockey"
        return "Class_I", "LINE", "Unknown"

    if any(x in name for x in ["hat", "tc1", "mariner", "piggybac", "mutator", "harbinger", "pif"]):
        if "hat" in name:
            return "Class_II", "TIR", "hAT"
        if "tc1" in name or "mariner" in name:
            return "Class_II", "TIR", "TcMar-Mariner"
        if "piggybac" in name:
            return "Class_II", "TIR", "PiggyBac"
        if "mutator" in name:
            return "Class_II", "TIR", "Mutator"
        if "harbinger" in name or "pif" in name:
            return "Class_II", "TIR", "PIF-Harbinger"
        return "Class_II", "TIR", "Unknown"

    if "sine" in name:
        return "Class_I", "SINE", "Unknown"

    return "Unknown", "Unknown", "Unknown"

def normalize_dfam_ct_terms(terms: list[str]) -> tuple[str, str, str]:
    terms_set = set(terms)

    dfam_class = "Unknown"
    order = "Unknown"
    superfamily = "Unknown"

    if "Class_I_Retrotransposition" in terms_set:
        dfam_class = "Class_I"

    elif "Class_II_Transposition" in terms_set:
        dfam_class = "Class_II"

    if "Long_Terminal_Repeat_Element" in terms_set:
        order = "LTR"

    elif "Non_LTR_Retrotransposon" in terms_set:
        order = "LINE"

    elif "Short_Interspersed_Element" in terms_set:
        order = "SINE"

    elif "Terminal_Inverted_Repeat_Element" in terms_set:
        order = "TIR"

    elif "Helitron" in terms_set:
        order = "Helitron"

    # Prefer the most specific final CT term as superfamily,
    # but avoid broad ontology terms.
    broad_terms = {
        "Interspersed_Repeat",
        "Transposable_Element",
        "Retrotransposon",
        "DNA_Transposon",
        "Class_I_Retrotransposition",
        "Class_II_Transposition",
        "Long_Terminal_Repeat_Element",
        "Non_LTR_Retrotransposon",
        "Short_Interspersed_Element",
        "Terminal_Inverted_Repeat_Element",
    }

    for term in reversed(terms):
        if term not in broad_terms:
            superfamily = term
            break

    # Normalize common Dfam labels.
    superfamily_map = {
        "Gypsy-ERV": "Gypsy",
        "Gypsy": "Gypsy",
        "Copia": "Copia",
        "BEL-Pao": "Pao",
        "Pao": "Pao",
        "CR1": "CR1",
        "RTE": "RTE",
        "R1": "R1",
        "R2": "R2",
        "Jockey": "Jockey",
        "L1": "L1",
        "Penelope": "Penelope",
        "hAT": "hAT",
        "TcMar": "TcMar",
        "PiggyBac": "PiggyBac",
        "Mutator": "Mutator",
        "PIF-Harbinger": "PIF-Harbinger",
        "Helitron": "Helitron",
    }

    superfamily = superfamily_map.get(superfamily, superfamily)

    return dfam_class, order, superfamily


def parse_dfam_hmm_metadata(hmm_path: str | Path) -> dict[str, dict]:
    """Read FamDB's nucleotide HMM metadata, including RM annotations and GA.

    Curated status follows Dfam accession namespaces (DF versus DR). Missing
    classification is retained as Unknown; model names are not taxonomy.
    """
    metadata = {}
    record = {}
    with open(hmm_path) as handle:
        for line in handle:
            fields = line.strip().split(maxsplit=1)
            if not fields:
                continue
            tag = fields[0]
            value = fields[1] if len(fields) > 1 else ""
            if tag in {"NAME", "ACC", "LENG", "ALPH", "GA", "TC", "NC", "CT"}:
                record[tag] = value
            elif tag == "CC":
                for label, key in (("Type:", "type"), ("SubType:", "subtype")):
                    if value.startswith(label):
                        record[key] = value[len(label):].strip()
            elif tag == "//":
                if "NAME" in record:
                    ct = record.get("CT", "")
                    cls, order, sf = normalize_dfam_ct_terms(
                        [term.strip() for term in ct.split(";") if term.strip()]
                    )
                    rm_type = record.get("type", "")
                    rm_map = {"LTR": ("Class_I", "LTR"), "LINE": ("Class_I", "LINE"),
                              "SINE": ("Class_I", "SINE"), "DNA": ("Class_II", "TIR"),
                              "RC": ("Class_II", "Helitron")}
                    if rm_type in rm_map:
                        cls, order = rm_map[rm_type]
                        sf = record.get("subtype") or "Unknown"
                    if sf == "Penelope":
                        cls, order = "Class_I", "Penelope"
                    accession = record.get("ACC", "")
                    meta = {"dfam_class": cls, "dfam_order": order,
                            "dfam_superfamily": sf, "dfam_ct": ct,
                            "accession": accession, "model_length": int(record.get("LENG", 0)),
                            "alphabet": record.get("ALPH", "").lower(),
                            "curated": True if accession.startswith("DF") else
                                       False if accession.startswith("DR") else None}
                    for key in ("GA", "TC", "NC"):
                        values = record.get(key, "").replace(";", "").split()
                        meta[key.lower()] = float(values[0]) if values else None
                    metadata[record["NAME"]] = meta
                record = {}
    return metadata
