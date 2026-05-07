# Heliconiini classification benchmark

Dataset: curated `Heliconius.lib`

| Evidence mode | Total compared | Class accuracy | Order accuracy | Superfamily accuracy | Notes |
|---|---:|---:|---:|---:|---|
| Pfam only | 538 | 0.3755 | 0.3662 | 0.1691 | Protein-domain evidence only |
| Dfam subset | 538 | 1.0000 | 1.0000 | 0.9405 | Includes Unknown expected labels |
| Dfam subset, ignore Unknown expected | 447 | 1.0000 | 1.0000 | 1.0000 | More informative benchmark |
| Dfam + heuristic TIR structure | 538 | 1.0000 | 1.0000 | 0.9405 | Same benchmark result; adds structure evidence table |
