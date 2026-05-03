# DRayTEannotator2

DRayTEannotator2 is a modular workflow for transposable element (TE) discovery, refinement, and genome annotation. It combines de novo repeat discovery, consensus extension, structural-aware curation, and post-annotation refinement into a reproducible and scalable pipeline.

The workflow builds upon the original DRayTE approach and extends it with modern Python-based orchestration, improved reproducibility, and post-processing modules inspired by RepeatCraft and EarlGrey.

## Workflow overview:

`Discovery → Extension → Reclassification → Curation → Annotation → Refinement → Summary`

Key outputs:

* Curated TE library (FASTA)
* Genome-wide TE annotation (RepeatMasker)
* Refined TE loci (GFF3/BED/TSV)
* Summary tables and publication-ready figures

⸻

Core pipeline

* RepeatModeler2-based de novo TE discovery
* RepeatMasker-based genome annotation
* Consensus extension (davidExtendConsRAM.pl)
* RepeatClassifier reclassification
* TE-Aid-assisted curation

New in DRayTEannotator2

1. Annotation refinement (TE-refine)

* RepeatCraft-inspired defragmentation
* Two modes:
    * strict: conservative merging
    * loose: permissive merging across intervening elements
* Nested TE detection
* Consensus-coordinate consistency filtering
* Outputs:
    * *.refinedRepeats.gff3
    * *.refinedRepeats.bed
    * *.annotation_refinement.tsv
    * *.annotation_refinement.stats.tsv
    * *.annotation_refinement.manifest.json

2. Nested-aware classification

Each locus is classified as:
```
LINE
LINE-nested
LTR
LTR-nested
...
```

This enables biologically meaningful separation of:

* primary insertions
* nested insertions

3. Summary generation (TE-summary)

* EarlGrey-like outputs
* High-level and family-level TE tables
* Divergence landscapes
* Superfamily profiles
* Publication-ready figures:
    * TE composition barplots
    * Repeat landscapes
    * Nested vs non-nested comparisons

⸻

## Installation

Recommended: Micromamba / Conda

Create environment:
```
micromamba create -n drayte -f environment.yml
micromamba activate drayte
pip install -e .
```

Alternative: Python virtualenv
```
python -m venv drayte
source drayte/bin/activate
pip install -e .
```

�

## Dependencies

Core

* Python ≥ 3.10
* pandas, numpy, matplotlib

TE discovery & annotation

* RepeatModeler2
* RepeatMasker
* RMBlast (recommended)
* RepeatClassifier
* RepeatScout

Additional tools

* MMseqs2
* DIAMOND
* CD-HIT
* EMBOSS (getorf)
* BLAST+
* BEDTools
* seqkit

Optional:

* HELIANO (Helitron detection)
* TE-Aid (manual curation)
* MMseqs2

### Required external scripts

The pipeline currently uses:

* `davidExtendConsRAM.pl` (included)
* RAMExtend
* associated RepeatModeler extension scripts

These should be available in the configured paths.

Recommended

* GNU Parallel
* Slurm (for cluster execution)
* `twoBitInfo` from KentUtils
* `faToTwoBit` from KentUtils

⸻

## Configuration

The pipeline runs using a YAML configuration file.

Example:

```
species: Cant
threads: 16
input_genome: genome.fa
output_dir: output/

repeatmodeler_dir: /path/to/RepeatModeler
repeatmasker_bin: RepeatMasker
repeatclassifier_bin: RepeatClassifier

# Extension
repeatmodeler_extend_script: davidExtendConsRAM.pl

# Refinement
refine_mode: loose
refine_max_gap: 150
```

Adjust all paths to your local installation.

⸻

## Running the Pipeline

Run with:

`drayte --config config.yaml`

⸻

## Post-processing

TE refinement
```
TE-refine \
  --rmout genome.out \
  --species Cant \
  --outdir Cant.annotation_refinement \
  --max-gap 150 \
  --mode loose
```

TE summary
```
TE-summary \
  --refined-tsv Cant.annotation_refinement/Cant.annotation_refinement.tsv \
  --genome genome.fa \
  --species Cant \
  --outdir Cant.summaryFiles
```

⸻


Output Structure


```
output/
├── discovery/
├── extension/
├── reclassify/
├── curation/
├── final_annotation/
├── annotation_refinement/
│   ├── *.gff3
│   ├── *.tsv
│   ├── *.stats.tsv
│   └── *.manifest.json
└── summaryFiles/
    ├── highLevelCount.txt
    ├── familyLevelCount.txt
    ├── divergence plots
    └── TE composition figures
```

## Important final output:

`curation/Final.RepeatModeler.Lib.fa`

This is the curated TE library for downstream RepeatMasker annotation.

⸻

## Recommended strategy

1. Run DRayTEannotator2 per species
2. Merge TE libraries across species
3. Run RepeatMasker with combined library
4. Apply TE-refine
5. Generate summaries with TE-summary

⸻

## Notes

Header format

RepeatMasker-compatible FASTA headers must preserve the canonical format:

`>family_id#CLASS/FAMILY`

Avoid non-standard formats such as:

`>family__CLASS___FAMILY#CLASS/FAMILY`

These can break downstream classification and summary parsing.

⸻

## Citation

If you use this workflow, please cite:

Novel sex-specific genes and diverse interspecific expression in the antennal transcriptomes of Ithomiine butterflies
F Cicconardi, BJ Morris, J Martelossi, DA Ray, SH Montgomery
Genome Biology and Evolution 16 (10), evae218
https://academic.oup.com/gbe/article-pdf/doi/10.1093/gbe/evae218/60605110/evae218.pdf

and ...
* RepeatModeler
* RepeatMasker
* TE-Aid
* original DRayTE workflow
* EarlGrey (if used in parallel)

and any additional tools used in your specific analysis.

⸻

## Contact

For issues, improvements, or collaboration:

GitHub Issues:

https://github.com/francicco/DRayTEannotator2
