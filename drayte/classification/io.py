import csv

from .models import Family


BOOL_FIELDS = {
    "ltr_present",
    "tir_present",
    "helitron_signal",
    "tsd_present",
    "polyA_present",
}

LEGACY_DOMAIN_BOOL_FIELDS = {
    "rt_present": "RT_domain",
    "integrase_present": "INT_domain",
    "transposase_present": "Transposase_domain",
}

INT_FIELDS = {
    "consensus_len",
    "n_copies",
    "orf_count",
    "orf_max_len",
}

FLOAT_FIELDS = {
    "homology_score",
    "boundary_consistency",
    "fragmentation_score",
}


def parse_bool(v):
    return str(v).strip().lower() in {"1", "true", "yes", "y"}


def parse_domains(v):
    if v in {None, "", "NA", "Unknown"}:
        return set()

    return {
        x.strip()
        for x in str(v).replace(",", ";").split(";")
        if x.strip()
    }


def load_families_tsv(path):
    families = []

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            clean = {}
            domains = parse_domains(row.get("domains"))

            for k, v in row.items():
                if k in LEGACY_DOMAIN_BOOL_FIELDS:
                    if parse_bool(v):
                        domains.add(LEGACY_DOMAIN_BOOL_FIELDS[k])
                    continue

                if k == "domains":
                    continue

                if k in BOOL_FIELDS:
                    clean[k] = parse_bool(v)

                elif k in INT_FIELDS:
                    clean[k] = int(v)

                elif k in FLOAT_FIELDS:
                    clean[k] = float(v)

                else:
                    clean[k] = v

            clean["domains"] = domains
            families.append(Family(**clean))

    return families


def write_classification_tsv(results, outpath):
    fields = [
        "family_id",
        "class",
        "order",
        "superfamily",
        "status",
        "confidence",
        "score",
        "evidence",
    ]

    with open(outpath, "w") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        for row in results:
            writer.writerow(row)
