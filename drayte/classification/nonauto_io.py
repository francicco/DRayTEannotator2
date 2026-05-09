from __future__ import annotations

import csv
from pathlib import Path

from .nonauto import NonAutonomousCall


def write_nonautonomous_tsv(
    calls: list[NonAutonomousCall],
    outfile: str | Path,
) -> None:

    fields = [
        "family_id",
        "candidate_type",
        "score",
        "confidence",
        "evidence",
    ]

    with open(outfile, "w") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()

        for call in sorted(
            calls,
            key=lambda x: (
                x.candidate_type,
                -x.score,
                x.family_id,
            ),
        ):
            writer.writerow({
                "family_id": call.family_id,
                "candidate_type": call.candidate_type,
                "score": call.score,
                "confidence": call.confidence,
                "evidence": call.evidence,
            })
