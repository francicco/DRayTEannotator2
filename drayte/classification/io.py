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
    "rt_present": "RT",
    "integrase_present": "INTEGRASE",
    "transposase_present": "TRANSPOSASE",
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


def write_evidence_tsv(families, results, outpath):
    result_by_id = {
        row["family_id"]: row
        for row in results
    }

    fields = [
        "family_id",
        "final_class",
        "final_order",
        "final_superfamily",
        "status",
        "confidence",
        "score",
        "homology_class",
        "homology_superfamily",
        "homology_score",
        "dfam_model",
        "dfam_order",
        "dfam_superfamily",
        "dfam_score",
        "domains",
        "ltr_present",
        "tir_present",
        "helitron_signal",
        "tsd_present",
        "polyA_present",
        "orf_count",
        "orf_max_len",
        "evidence",
    ]

    with open(outpath, "w") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        for f in families:
            result = result_by_id[f.family_id]

            writer.writerow({
                "family_id": f.family_id,
                "final_class": result["class"],
                "final_order": result["order"],
                "final_superfamily": result["superfamily"],
                "status": result["status"],
                "confidence": result["confidence"],
                "score": result["score"],
                "homology_class": f.homology_class,
                "homology_superfamily": f.homology_superfamily,
                "homology_score": f.homology_score,
                "dfam_model": f.dfam_model or "",
                "dfam_order": f.dfam_order or "",
                "dfam_superfamily": f.dfam_superfamily or "",
                "dfam_score": f.dfam_score,
                "domains": ";".join(sorted(f.domains)),
                "ltr_present": f.ltr_present,
                "tir_present": f.tir_present,
                "helitron_signal": f.helitron_signal,
                "tsd_present": f.tsd_present,
                "polyA_present": f.polyA_present,
                "orf_count": f.orf_count,
                "orf_max_len": f.orf_max_len,
                "evidence": result["evidence"],
            })
