# Optional final Dfam validation

The nested `classification.dfam` configuration adds conservative reference-HMM
validation after structural, protein-domain and MMseqs classification evidence
has been collected. It reuses `classification/dfam.py` for NHMMER searches and
`classify_family` for confidence scoring. It does not add Dfam sequences to the
species-specific library or change the final RepeatMasker command.

## Configuration

Merge [the example](../config/dfam-validation.example.yaml) into your pipeline
YAML. The database may be a FamDB directory or an exported nucleotide HMM file.
There are no built-in database paths or automatic reference downloads.

```yaml
classification:
  dfam:
    enabled: true
    database: /path/to/Dfam_4.0/FamDB
    taxon: Arthropoda
    curated_only: true
    engine: nhmmer
    threads: 16
```

`taxon` can be an unambiguous name such as Arthropoda, Lepidoptera or Collembola,
or an NCBI taxonomy identifier. Ancestors and descendants are included by default;
either can be disabled with `ancestors: false` or `descendants: false`. Necessary
FamDB partitions must be installed. No hardcoded partition numbers are used.

For an existing HMM export, use its path and `taxon: null`. That file is treated
as already taxonomically selected; passing a taxon with an HMM file raises an
error rather than silently pretending to filter it. Curated-only acceptance still
applies to exported HMM hits using Dfam's DF/DR accession namespaces. Set
`curated_only: false` to also consider uncurated models when available.

When this nested mapping is present, it supersedes the legacy flat `dfam_db`
configuration and existing legacy Dfam result caches. There is one Dfam search.
`enabled: false` bypasses Dfam entirely, including any old Dfam cache. Omitting
the nested mapping retains the existing legacy integration behavior.

## Dependencies and FamDB export

Install HMMER (including `nhmmer`) and the official FamDB CLI matching your
reference format. `famdb_bin` and `nhmmer_bin` accept executable paths. FamDB's
Python dependencies are available through:

```sh
pip install -e '.[famdb]'
```

This extra installs h5py and NumPy; the FamDB CLI itself and HMMER must be installed
separately. They are unnecessary when validation is disabled or, for FamDB only,
when using an existing HMM export.

Before export, DRayTE checks the installed CLI help and captures database `info`.
The export uses the verified CLI shape:

```sh
famdb.py -i /path/to/FamDB families --format hmm \
  --ancestors --descendants --curated Arthropoda
```

Official references:
[FamDB](https://github.com/Dfam-consortium/FamDB),
[Dfam 4.0 partitions](https://dfam.org/releases/Dfam_4.0/families/FamDB/), and
[HMMER](http://hmmer.org/).

FamDB owns the HDF5 layout; DRayTE does not read partitions directly. Empty
exports and export diagnostics are treated as errors so missing partitions or
unresolved taxa cannot silently produce a successful empty validation. FamDB
versions with different CLI capabilities fail with an actionable message.
The parser reads `ACC`, `LENG`, `ALPH`, `GA`, `TC`, `NC`, `CT`, and the exported
RepeatMasker Type/SubType annotations. Nucleotide HMMs are checked explicitly;
protein HMMs for `hmmscan` are rejected. Family-name guesses are not used in the
new validation decisions.

`on_missing: error` is the default. `on_missing: warn` explicitly permits unavailable
references/tools, failed CLI commands or empty selections to continue with the
baseline classification and a recorded `unavailable` validation status. Invalid
configuration or incompatible model metadata still raises an error. Unavailable
runs are retried on resume, rather than recorded as successful validation.

## Decisions and confidence

The existing classifier first produces a baseline without Dfam evidence. Reliable
Dfam hits must pass all of:

- E-value <= `max_evalue` (default 1e-5) and score >= `min_score` (20).
- The model's GA gathering threshold, when provided. With `require_ga: true`
  (default), models without GA cannot affect classification. TC and NC are
  recorded as metadata, not substituted for GA. FamDB's ordinary `hmm` export
  supplies general cutoffs; species-specific TH thresholds are not applied to a
  broad clade. Disabling `require_ga` explicitly permits the global-score fallback.
- Aligned consensus coverage >= `min_query_coverage` (0.10), and model coverage
  >= `min_model_coverage` (0.50). Coverage uses inclusive coordinates and handles
  either strand. Consensus coverage is the fraction of the input TE consensus;
  NHMMER itself calls the HMM its query. Each alignment is evaluated individually,
  without stitching separate fragments into apparent full-length evidence.
- Mapped class and order, and curated status when requested.

These are conservative configurable policy defaults, not empirically calibrated
accuracy estimates. Short conserved motifs alone should not determine a family.

Matches scoring at least `competition_score_fraction` (0.95) times the best
accepted hit compete with it. Contradictory known class, order or superfamily
labels make the assignment ambiguous and block transfer. This score comparison
is a competition heuristic, not a probability of correctness. Common aliases
such as Bel-Pao/Pao and PiggyBac/piggyBac are compared consistently.

Reliable unambiguous evidence can fill unresolved hierarchy nodes. Dfam-only
assignments remain `putative`, with at most MEDIUM confidence and score 0.64.
Any conflicting existing known node is retained and reported, including a
supported structural or protein-domain assignment. A weak or absent hit leaves
Unknown unchanged. Existing superfamily assignments are not replaced.

`independent: false` is the default: concordant Dfam evidence is recorded without
an automatic confidence increase, because RepeatClassifier headers and nucleotide
reference FASTAs can already derive from Dfam. Only set `independent: true` after
establishing that the baseline homology is independent. This explicit assertion
allows the existing classifier's score/confidence result to support a concordant
assignment; it cannot override conflicts, and Dfam-only calls remain provisional.
The pipeline does not infer independence from different search algorithms.

## Reports and final library

Existing classification and evidence TSVs retain their usual columns and add
original class/order/superfamily, confidence, score and supporting evidence,
plus `dfam_agreement` and `dfam_contributed`, when validation is enabled.

Additional files follow the species-prefix convention:

- `<species>.dfam_evidence.tsv`: every reported alignment, including rejected
  matches, accession/name, mapped hierarchy, original CT, curation status,
  score/E-value, model threshold, lengths/coordinates, both coverage measures,
  rejection reason, agreement, contribution and original/final labels/confidence.
- `<species>.dfam_conflicts.tsv`: conflicting or competing ambiguous accepted
  matches. An empty conflict table still has a header.
- `dfam_validation/dfam.reference.json`: source, FamDB database information and
  export command when applicable, exact exported HMM digest and model count.
- `dfam_validation/dfam.validation.manifest.json`: effective policy, tool help
  and executable digest, input digests, search size and result location.

The final `<species>.classified.fa` is written by the existing RepeatMasker
header formatter, retaining original family identifiers and consensus bases.
It uses established labels such as `family#LINE/CR1` and `family#DNA/TcMar-Mariner`.
The existing final annotation still uses this library with `RepeatMasker -lib`;
there is no concatenation of Dfam consensuses and no Dfam-only genome annotation.

## Restart behavior

Selection/export, raw searches and classification policy have separate provenance.
Changing taxonomy or curation selection refreshes the export; changed exported
model bytes or NHMMER executable select a new raw-search cache. Identical model
content can reuse the search even when selection settings differ. Changing coverage, acceptance or competition
thresholds reuses compatible NHMMER chunks and regenerates classification reports.
Changing Dfam configuration does not force discovery or consensus extension.
Other existing biological-input mismatch checks remain active.

NHMMER uses `-Z` with the full input-consensus database size, so changing chunk
boundaries does not change the search-space size used for E-values. Completed
chunks require a matching manifest and NHMMER completion marker. Partial chunks
are rerun. Compatible historical chunk partitioning is retained on resume.

FamDB partition names/sizes/modification times plus database information and CLI
identity invalidate exports; the exported HMM content is hashed. Preserve the
source partitions and provenance for reproducibility. In-place partition edits
that deliberately preserve both file size and modification time require removing
the export manifest to force a fresh export. Existing HMM exports are hashed on
every validation attempt. Input files changing during execution prevent a
completion from being recorded.

## Standalone classifier

Put only the contents of the `dfam` mapping in `dfam-options.yaml`, then run:

```sh
drayte-classify --fasta TE_library.fa --pfam-db Pfam-A.hmm \
  --dfam-config dfam-options.yaml --dfam-outdir classification_dfam \
  --output classifications.tsv --evidence-output evidence.tsv
```

This mode requires FASTA input and cannot be combined with legacy `--dfam-db` or
`--dfam-tblout`. Standalone Dfam report prefixes use the output TSV stem.

## Validation performed

Unit and pipeline tests cover parsing, hierarchy mapping, acceptance/coverage,
weak and competing hits, conflict protection, disabled behavior, missing tools
and references, output format and policy/database changes on resume.

A real compatibility test was also run with HMMER 3.4, official FamDB 3.0.0
(source revision `630f0d7881d021e70d7fcfeaac820eaa79a787fb`), and a synthetic
FamDB HDF5 database holding the historical MADE1 model distributed with HMMER.
Its added GA is a test cutoff. This exercises actual export, NHMMER parsing,
provisional assignment, a negative control and completed-chunk reuse; it is not
a biological validation against the full Dfam 4.0 release.

To rerun that optional test against an installed compatible fixture/reference:

```sh
DRAYTE_FAMDB_TEST_DATABASE=/path/to/test/FamDB \
DRAYTE_FAMDB_BIN=/path/to/famdb.py \
DRAYTE_NHMMER_BIN=/path/to/nhmmer \
python -m pytest -q tests/test_dfam_external_tools.py
```

The reference must include a classified MADE1 nucleotide model with GA; taxon 1
is used unless `DRAYTE_FAMDB_TEST_TAXON` is set. Ordinary test runs skip this
external test when its tools/database are not configured.
