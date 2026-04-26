# DRayTEannotator2

DRayTEannotator2 is a transposable element (TE) annotation workflow designed for de novo repeat discovery, consensus extension, reclassification, curation, and final TE library generation for downstream genome annotation.

The pipeline is inspired by the original DRayTE workflow and modernized into a modular Python-based implementation with clearer stage separation, reproducibility, and easier maintenance.

The standard workflow is:

`Discovery → Extension → Reclassify → Curation → Report`

The final curated TE library can then be used for genome-wide annotation with RepeatMasker or integrated with TE libraries from related species for iterative masking and comparative analyses.

⸻

## Main Features

* RepeatModeler-based de novo TE discovery
* RepeatMasker-based first-pass annotation
* Consensus extension using davidExtendConsRAM.pl
* RepeatClassifier reclassification
* TE-Aid assisted curation support
* Final curated TE library generation
* Designed for HPC / cluster environments
* Config-driven execution
* Modular Python implementation replacing older shell-heavy workflows

⸻

## Dependencies

Core software

The following tools must be installed and accessible.

Required

* Python >= 3.9
* RepeatModeler 2
* RepeatMasker
* RMBlast (recommended) or CrossMatch
* RepeatClassifier
* RepeatScout
* Perl
* BioPerl
* BEDTools
* BLAST+
* seqtk
* samtools
* cd-hit
* EMBOSS
* MAFFT
* TE-Aid

Required external scripts

The pipeline currently uses:

* davidExtendConsRAM.pl
* RAMExtend
* associated RepeatModeler extension scripts

These should be available in the configured paths.

Recommended

* GNU Parallel
* Slurm (for cluster execution)
* twoBitInfo from KentUtils
* faToTwoBit from KentUtils

⸻

## Python Environment

Create a dedicated environment:

python -m venv DRayTEannotator2
source DRayTEannotator2/bin/activate
pip install -e .

If using Conda or Micromamba:

micromamba create -n DRayTEannotator2 python=3.10
micromamba activate DRayTEannotator2
pip install -e .

⸻

## Installation

Clone the repository:

```git clone https://github.com/francicco/DRayTEannotator2.git
cd DRayTEannotator2
```

Install the package:

`pip install -e .`

This will make the drayte command available.

⸻

## Configuration

The pipeline runs using a YAML configuration file.

Example:

species: Cant
threads: 16
input_genome: /path/to/genome.fasta
output_dir: /path/to/output
repeatmodeler_dir: /path/to/RepeatModeler
repeatmasker_bin: /path/to/RepeatMasker
repeatclassifier_bin: /path/to/RepeatClassifier
repeatscout_dir: /path/to/RepeatScout
repeatmodeler_extend_script: /path/to/davidExtendConsRAM.pl
te_aid_dir: /path/to/TE-Aid
blast_bin: /path/to/blast+

Adjust all paths to your local installation.

⸻

## Running the Pipeline

Run with:

`drayte --config config.yaml`

⸻

Output Structure

Typical output layout:

```output/
├── discovery/
│   ├── assemblies_dir/
│   ├── rmodeler_dir/
│   └── rmasker_dir/
│
├── extensions/
│   ├── extract_align/
│   ├── extensionwork/
│   └── extendlogs/
│
├── reclassify/
│
├── curation/
│   ├── prioritize/
│   ├── te-aid/
│   └── Final.RepeatModeler.Lib.fa
│
└── reports/
```

## Important final output:

curation/Final.RepeatModeler.Lib.fa

This is the curated TE library for downstream RepeatMasker annotation.

⸻

## Recommended Final Annotation Strategy

For best results:

1. Run DRayTEannotator2 independently for each species
2. Collect curated TE libraries
3. Merge libraries across related species
4. Remove redundancy if necessary
5. Perform final RepeatMasker annotation using the combined library

This improves sensitivity and captures lineage-specific TE diversity.

⸻

## Notes

Header format

RepeatMasker-compatible FASTA headers must preserve the canonical format:

`>family_id#CLASS/FAMILY`

Avoid non-standard formats such as:

`>family__CLASS___FAMILY#CLASS/FAMILY`

These can break downstream classification and summary parsing.

⸻

## EarlGrey integration

EarlGrey can be run in parallel as an independent validation workflow, particularly useful for:

* Helitron detection
* TIR detection
* LTR structural validation
* comparison of conservative vs permissive TE discovery

Current recommendation:

Use DRayTE as the main discovery workflow and EarlGrey as complementary structural validation.

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
