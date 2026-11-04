# DRayTEannotator2

DRayTEannotator2 is a modular framework for transposable element (TE) discovery, consensus reconstruction, evidence-based classification, genome annotation, and post-annotation refinement.

The workflow formalises and extends the original DRayTE strategy developed in the Ray laboratory into a reproducible and scalable Python-based framework. DRayTEannotator2 integrates de novo repeat discovery, consensus extension, structural evidence, protein-domain annotation, reference-assisted similarity rescue, and nested-aware post-processing into a unified workflow for TE annotation in eukaryotic genomes.

The pipeline is designed both for:
- high-quality TE library reconstruction
- genome-wide TE annotation and quantification

and is particularly suitable for non-model genomes where repeat families are fragmented, poorly classified, or structurally under-represented.

---

## Workflow overview

```text
Discovery
    ↓
Consensus extension
    ↓
Reclassification
    ↓
Evidence-based curation
    ↓
Optional structure-aware discovery
    ↓
Genome annotation
    ↓
TE-refine post-processing
    ↓
Summary generation

-----

## TE discovery and reconstruction

* RepeatModeler2-based de novo repeat discovery
* RepeatMasker-compatible consensus normalisation
* Consensus extension using:
    * Extract_ALIGN
    * davidExtendConsRAM.pl
    * RAMExtend
* Recovery of fragmented or truncated TE consensuses
* Multi-copy consensus reconstruction from genomic context

### New in DRayTEannotator2

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

----- 

## TE-refine post-processing

DRayTEannotator2 implements overlap-aware TE refinement inspired by RepeatCraft and interval-graph optimisation strategies.

Features include:

* RepeatMasker defragmentation
* nested insertion detection
* overlap-aware annotation
* conflict resolution
* interval scheduling optimisation
* divergence-aware scoring

Two complementary annotation layers are produced:

* overlap-aware annotations
* filtered non-overlapping annotations

Summary generation

The pipeline generates publication-ready TE summaries including:

* class-level TE composition
* family-level abundance
* divergence landscapes
* nested vs non-nested TE profiles
* genome coverage statistics
* repeat landscapes

Outputs are conceptually comparable to EarlGrey summary statistics while retaining DRayTE-specific refinement logic.

-----

## Evidence-based TE classification

DRayTEannotator2 implements a multi-evidence TE classification framework integrating:

* RepeatClassifier annotations
* Dfam nhmmer searches and optional [final Dfam/FamDB validation](docs/dfam-validation.md)
* Pfam protein-domain annotation
* ORF prediction
* structural signatures
* MMseqs2 reference rescue using one-way similarity
* RepeatMasker header homology
* genomic copy support

1. Classification is hierarchical and confidence-aware.

Assignments are evaluated at:

* TE class
* TE order
* TE superfamily

2. Confidence levels are calibrated according to agreement among independent evidence sources.

3. Reference-assisted `MMseqs2` rescue

4. Unclassified families can be rescued using nucleotide similarity searches against classified TE libraries; one-way transfer allows several predicted families to match the same reference.

5. Rescue filtering includes:

* minimum sequence identity
* minimum alignment length
* minimum query coverage
* minimum target coverage
* reference-label agreement filtering
* self-hit exclusion

6. Hierarchical rescue is supported:

* order-level rescue without forcing uncertain superfamily assignments

7. Structure-aware discovery

Optional structure-based TE detection includes:

* HELIANO integration for Helitron discovery

8. Structurally detected candidates are:

* normalised to RepeatMasker-compatible format
* compared against curated libraries using MMseqs2
* dereplicated before integration

This enables recovery of low-copy or structurally defined elements not captured efficiently by repeat abundance alone.

-----

## Scalable execution

Large TE libraries can be processed using chunked parallel execution for:

* hmmscan (Pfam)
* nhmmer (Dfam)

Features include:

* automatic chunk splitting
* parallel execution
* cached intermediate reuse
* automatic regeneration of missing files
* restart-safe workflows

Intermediate outputs are reused automatically across reruns.

-----

For the pinned RepeatMasker 4.2.4 / RepeatModeler 2.0.9 / FamDB 3.0.0
cluster environment, Dfam 4.0 data installation, external `cross_match`
requirement and configuration paths, see [the environment guide](docs/environment.md).
Installing the software does not download the Dfam database.

## Installation

Recommended: Micromamba / Conda

Create environment:
```
micromamba create -n drayte -f environment.yml --strict-channel-priority
micromamba activate drayte
pip install -e .
```

Alternative: Python virtualenv
```
python -m venv drayte
source drayte/bin/activate
pip install -e .
```

-----

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

-----

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

-----

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
extension_workers: null  # defaults to threads; use 1 for serial family extension
repeatmodeler_extend_script: davidExtendConsRAM.pl
ucsc_tools_dir: /path/to/kent/bin
repeat_after_me_dir: /path/to/RepeatAfterMe
phrap_dir: /path/to/phrap

# Refinement
refine_mode: loose
refine_max_gap: 150
```

Adjust all paths to your local installation.

Consensus extension runs independent TE families concurrently, with one external
extension process per worker. `extension_workers` must be a positive integer, or
`null` (the default) to use `threads`. For example, `threads: 16` with
`extension_workers: 4` uses 16 BLASTN threads in the preceding search and up to four
concurrent family-extension jobs afterward. Use `extension_workers: 1` for the
previous serial behavior. The worker count is capped at the number of families.

Each family retains its own working directory and log under
`extension/extensionwork/` and `extension/extendlogs/`. Combined extension log
messages include the family name, and the summary keeps extraction order even
when jobs finish out of order. Individual failures are recorded without stopping
other families. The external extension commands and biological settings are
unchanged; workers control family concurrency, not threading inside those tools.
Reduce workers if memory or disk access becomes limiting, and allow for any
additional CPUs used internally by your installed tools.

Changing only `extension_workers` preserves existing checkpoints: rerun the same
pipeline command to resume, skipping completed families and retrying failed or
unfinished families using the existing resume checks. A completed extension stage
is reused without rerunning it.

-----

## Running the Pipeline

Run with:

`drayte --config config.yaml`

Supply a reference FASTA for classification rescue on the command line with:

```bash
drayte --config config.yaml \
  --mmseqs-rescue-db /path/to/Heliconius.classification.fa
```

This flag overrides `classification.mmseqs_rescue_db` in the YAML file while
preserving the other classification settings. Without the flag, the YAML value
is used. It selects classification evidence; it does not add the reference
sequences directly to the final annotation library. The same option is also
available on the standalone `drayte-classify` command.

External-reference rescue uses one-way query-to-reference search. Several
predicted consensuses can receive classification evidence from the same curated
reference. Identity, alignment length and coverage filters still apply, and
competitive reference matches must agree on the transferred classification.
The near-best score window defaults to 95% of the best bit score and can be changed
with `--mmseqs-rescue-score-fraction`.
Query and reference identifiers remain separate, and reference sequences are
never added to the output family set. The search cache fingerprints both FASTAs
by content, parameters, policy and algorithm version. Logs distinguish alignments,
families with rescue evidence and labels actually changed by that evidence. See
[Reference-assisted classification](docs/mmseqs-reference-rescue.md) for outputs,
thresholds, interpretation and resuming an older assisted run.

Rerunning the same command against the same output directory automatically
resumes the analysis. A stage is skipped only after its manifest and required
outputs have been validated. For extension, the cached query is retained as the
downstream structure-analysis base so a historical checkpoint is not mixed with
a subsequently regenerated discovery library.

Classification checks the resolved structural-input content rather than treating
the internal `structure_base_fasta` path added during extension resume as a new
biological setting. The derived structural-input FASTA is reconstructed from the
selected sources before fingerprinting, so older source timestamps cannot hide
changed content. Older classification checkpoints are migrated only after
their input fingerprints, applicable settings hashes, and required outputs have
been validated. Changed inputs or classification settings still stop resume;
the error now names the differing fingerprint or setting. The optional nested
`classification.dfam` validation settings are refreshed separately: changing
Dfam selection, reference contents, or validation thresholds rebuilds the affected
search/decisions while retaining compatible upstream evidence.

Use a new output directory when intentionally starting an analysis with different
inputs or settings.

RepeatClassifier can legitimately transform consensus bases while preserving
family IDs and lengths. Resume validates its completion checkpoint rather than
requiring output sequences to equal the unclassified input. Legacy results can
be migrated with a successful stage receipt and matching hashes, or consistent
input/output/receipt timestamps when historical hashes are absent. Identical
consensus merges preserve the existing merged file and timestamp. Invalid modern
checkpoints never fall back to legacy adoption, and a new attempt is recorded as
`running` until successful validation. Discovery also skips database rebuilding
when a usable cached RepeatModeler library is already available.

If older retries already invalidated the legacy timestamps, a completed curated
library and its family table may support one-time recovery. Preview with
`python -m drayte.pipeline.reclassify_recovery --config config.yaml`; this does
not run RepeatClassifier or modify output files. See
[RepeatClassifier recovery](docs/repeatclassifier-recovery.md) for validation,
backups, and the apply/resume procedure.

### Storage management for many genomes

Use `drayte --config config.yaml --storage-mode compact` to remove validated
completed-stage scratch and archive bulky evidence. Completed extension scratch
is compressed before classification, retaining the consensus, alignment, image,
and REPAM files needed by downstream stages. Remaining evidence is archived after
the whole run succeeds; unfinished extension work stays loose.
The default remains `keep`. The equivalent YAML is `storage: {mode: compact}`.
Core final libraries, RepeatMasker `.out/.tbl`, refined annotations, summaries,
and completion manifests remain directly accessible. Per-family work, selected
search chunks, masked FASTAs, and final alignments are stored in verified
compressed archives under `.storage/`; resuming `drayte` restores them before
checking stage fingerprints. Do not delete that directory.

For an existing **idle** run, preview with:

```bash
python -m drayte.storage.cli --config Eisa.config.yaml
```

Add `--apply` to compact, or `--restore --apply` to restore loose evidence for
inspection. After reinstalling the editable package, `drayte-storage` is an
equivalent command. Full procedures, retained files, archive verification,
locking, temporary quota headroom, and legacy checkpoint limits are documented
in [Storage management](docs/storage-management.md).

RepeatMasker's full stdout and stderr are saved to files rather than printed in
the terminal or copied into `drayte.log`:

- Initial pass: `discovery/rmasker_dir/<species>.RepeatMasker.log`
- Final pass: `final_annotation/<species>.RepeatMasker.log`

The main log shows brief start, completion, or failure messages with the full-log
path. Each rerun appends an attempt header with a UTC timestamp and command,
preserving earlier output. Cache skips leave existing logs untouched, and older
completed runs without these logs remain valid. Available log paths are also
recorded as `repeatmasker_log` in the corresponding stage manifest.

### Standalone EarlGrey storage utility

[`compact_earlgrey.py`](compact_earlgrey.py) is a single file that can be copied
anywhere and run with standard-library Python 3.10+. It needs no DrayAnn
installation, `pip` packages, or Conda environment. It preserves final summaries,
libraries, and raw RepeatMasker `.out/.tbl` while reversibly archiving bulky work:

```bash
python3 /path/to/compact_earlgrey.py --root Avcr_Earlgrey
# Only after that output directory is idle and the preview is eligible:
python3 /path/to/compact_earlgrey.py --root Avcr_Earlgrey --apply --idle
```

See [EarlGrey storage](docs/earlgrey-storage.md) for batch previews, restoration,
and wrapping future EarlGrey commands to restore before resuming and compact
after successful completion. Native EarlGrey calls require an explicit restore
first. The script also reads archives created by the previous packaged utility.

-----

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

-----

## Output Structure


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

-----

## Important final output:

`curation/Final.RepeatModeler.Lib.fa`

This is the curated TE library for downstream RepeatMasker annotation.

-----

## Recommended strategy

1. Run DRayTEannotator2 per species
2. Merge TE libraries across species
3. Run RepeatMasker with combined library
4. Apply TE-refine
5. Generate summaries with TE-summary

-----

## TE classification

Evidence-based classifier

Example:

```
drayte-classify \
  --fasta TE_library.fa \
  --pfam-db Pfam-A.hmm \
  --dfam-db dfam.hmm \
  --pfam-outdir classification_pfam \
  --dfam-outdir classification_dfam \
  --mmseqs-rescue-db TE_library.fa \
  --chunks 40 \
  --jobs 40 \
  --cpu 1 \
  --output classifications.tsv \
  --evidence-output evidence.tsv
```

For conservative final validation from Dfam 4.0 FamDB, use the nested
`classification.dfam` configuration in
[config/dfam-validation.example.yaml](config/dfam-validation.example.yaml).
It supersedes the legacy `dfam_db` search, records original and final assignments,
and reports conflicts. Dfam reference sequences are never appended to the final
species-specific RepeatMasker library. See the
[validation guide](docs/dfam-validation.md) for taxonomy, thresholds, dependencies,
standalone CLI usage, and resume behavior.

Outputs:

* classification table
* evidence table
* domain annotations
* reference-assisted rescue evidence

-----

# Preparing Pfam and Dfam databases

DRayTEannotator2 supports:
- Pfam protein-domain annotation via `hmmscan`
- Dfam TE homology annotation via `nhmmer`

Both databases must be prepared before running `drayte-classify`.

---

# Pfam setup

## Download Pfam

Download the current Pfam-A HMM database from the official Pfam/InterPro FTP server:

```bash id="v4p1js"
mkdir -p Databases/PFAM
cd Databases/PFAM

wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
```

---

## Prepare Pfam database

Press the HMM database for efficient `hmmscan` execution:

```bash id="q1yn0w"
hmmpress Pfam-A.hmm
```

This generates:

```text id="j3hnm4"
Pfam-A.hmm.h3f
Pfam-A.hmm.h3i
Pfam-A.hmm.h3m
Pfam-A.hmm.h3p
```

Required tool:
- HMMER3

---

## Example usage

```bash id="vbb2z3"
drayte-classify \
  --fasta TE_library.fa \
  --pfam-db /path/to/Pfam-A.hmm \
  --pfam-outdir classification_pfam \
  --output classifications.tsv
```

---

# Dfam setup

## Configure Dfam 4.0 validation

[DEFAULT.config.yaml](DEFAULT.config.yaml) contains the active nested
`classification.dfam` block. Set its `database` to your installed FamDB directory
and its `taxon` to the study group. The template uses Collembola, matching its
preset; use Lepidoptera for a Heliconius run. Set `enabled: false` to bypass Dfam.

The pipeline exports selected nucleotide profiles with the installed FamDB CLI.
For a manual export, the equivalent command is:

```bash
famdb.py -i /path/to/FamDB families --format hmm \
  --ancestors --descendants --curated Arthropoda > dfam-arthropoda.hmm
```

With an existing exported HMM file, set `database` to that file and `taxon: null`.
Renaming a full Dfam library does not select an arthropod subset. NHMMER reads
this ASCII HMM library directly; `hmmpress` is not required for this search.

See [the validation guide](docs/dfam-validation.md) for dependencies, thresholds,
reports and a standalone cluster test using an existing consensus FASTA. Dfam
reference sequences are not added to the final species-specific masking library.

---

# Recommended strategy

For most analyses:

- Pfam is used to detect conserved TE protein domains
- Dfam is used for nucleotide-level TE homology
- MMseqs2 rescue is used for recovering unclassified families through reference similarity

Using all three evidence layers substantially improves:
- TE order assignment
- superfamily resolution
- recovery of fragmented consensuses
- classification confidence

---

# Parallel execution

Both Pfam and Dfam searches support chunked parallel execution:

```bash id="x4b5vn"
--chunks 40 --jobs 40 --cpu 1
```

Where:
- `chunks` = number of FASTA partitions
- `jobs` = maximum concurrent jobs
- `cpu` = CPUs per individual hmmscan/nhmmer process

This is particularly important for large TE libraries.

---

# Cached execution

DRayTEannotator2 automatically reuses existing outputs:

Pfam cache:
```text id="vl9f6s"
classification_pfam/domains.domtblout
```

Dfam cache:
```text id="aq6d7m"
classification_dfam/dfam.merged.tblout
```

If cached outputs exist, searches are skipped automatically.

-----

## Downloading lineage-specific Dfam databases

DRayTEannotator2 supports lineage-specific Dfam libraries for more sensitive TE homology detection.

Users can generate custom Dfam HMM libraries directly from the Dfam browser.

### Step 1 — Open the Dfam browser

Go to:

[Dfam Browse Page](https://dfam.org/browse)

### Step 2 — Select the taxonomic group

In the `Taxon` field:

- type the lineage of interest

- select the desired taxon

- optionally enable `Descendants`

Example for Arthropoda:

[Arthropoda Dfam Browser Example](https://dfam.org/browse?clade=6656&clade_descendants=true)

This retrieves all Dfam entries associated with:

- Arthropoda

- descendant clades

This is generally recommended for TE annotation in non-model arthropods.

### Step 3 — Download the HMM library

At the bottom of the page, click:

```text id="t4cf2y"

HMM

```

This downloads the lineage-specific Dfam HMM database.

### Step 4 — Prepare the database

After download:

```bash id="d1f8za"

gunzip downloaded_file.hmm.gz

hmmpress downloaded_file.hmm

```

This generates the HMMER index files required for `nhmmer`.

-----

## Comparing EarlGrey and DRayAnn annotations

`compare_earlgrey_drayann.py` discovers matching species directories, reads the
two pipelines' `highLevelCount.txt` files, and writes paired tables and plots.
For example:

```bash
python compare_earlgrey_drayann.py \
  --root /user/work/tk19812/HoneyPotProject/TEannotation \
  --outdir /user/work/tk19812/HoneyPotProject/TEannotation/EarlGrey_vs_DRayAnn
```

The default uses the DRayAnn filtered annotation. Use
`--drayann-summary unfiltered` or `--drayann-summary both` to inspect the other
annotation layer. Existing nested rows are preserved in the raw output and
also collapsed into their parent TE classes for paired comparisons.

The comparison is descriptive rather than an accuracy benchmark: without an
independent truth set, greater coverage or copy number does not necessarily
mean a more accurate annotation.

-----

## Relationship to EarlGrey and RepeatCraft

DRayTEannotator2 is complementary to EarlGrey and RepeatCraft.

EarlGrey

EarlGrey provides:

* integrated TE annotation
* standardised reporting
* repeat landscapes
* automated summaries

DRayTEannotator2 emphasises:

* consensus reconstruction
* evidence-aware classification
* library refinement
* reference rescue using one-way similarity
* structural integration
* post-annotation optimisation

RepeatCraft

RepeatCraft focuses primarily on:

* RepeatMasker defragmentation
* nested repeat reconstruction

DRayTEannotator2 extends this concept with:

* interval-graph optimisation
* divergence-aware conflict resolution
* overlap-aware vs filtered annotation layers
* integrated refinement statistics

-----

## Notes

Header format

RepeatMasker-compatible FASTA headers must preserve the canonical format:

`>family_id#CLASS/FAMILY`

Avoid non-standard formats such as:

`>family__CLASS___FAMILY#CLASS/FAMILY`

These can break downstream classification and summary parsing.

-----

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

-----

## Contact

For issues, improvements, or collaboration:

GitHub Issues:

https://github.com/francicco/DRayTEannotator2
