# Reference-assisted classification

Use a classified TE FASTA as homology evidence for the families being classified:

```bash
drayte --config Hmel.drayannConfig.Assisted.yaml \
  --repeatmodeler-library "$RModelerFasta" \
  --mmseqs-rescue-db "${HELICONIINILIB%.lib}.classification.fa"
```

The CLI flag overrides `classification.mmseqs_rescue_db` in the configuration.
The reference headers need RepeatMasker class labels, for example
`>Helitron-1_Hmel#RC/Helitron`. The reference supplies classification evidence;
its sequences are not appended to the final annotation library. The existing
protein-domain, Dfam, structural and header evidence still participates in the
classification rules. Rescue does not guarantee a changed label for every hit,
and a changed label is not independently established biological accuracy.

## Reciprocal searches

The external-reference path runs two nucleotide searches:

1. Classification-query FASTA against the reference.
2. Reference against that same query FASTA.

A forward match must be the best retained match for its query, and its reference
target must have that query as its best reverse match. Query and reference
namespaces remain distinct, including when the same family name occurs in both
libraries. Reverse hits supply support for forward matches; they cannot create
new reference-only families in the classification output.

The search prefilter remains at 75% identity and 50% coverage. Existing retained
alignment filters remain at 80% identity, 500 aligned bases, 40% query coverage
and 30% reference coverage. The reverse search applies the corresponding
coverage thresholds with the library roles swapped. Assigning a rescued
superfamily still requires at least 80% identity and 80% coverage of both
sequences. The integrated pipeline can override retained-alignment thresholds
with `classification.mmseqs_min_identity`, `mmseqs_min_aln_len`,
`mmseqs_min_query_cov` and `mmseqs_min_target_cov`.

Strict reciprocal matching is conservative: several predicted families can
match one reference family, while only one is its best reverse match. Conflicting
or stronger existing evidence can also retain the original classification.
Report reference coverage and classification agreement separately when
benchmarking. Using the same reference as both pipeline input and evaluation
reference does not provide an independent accuracy estimate.

## Search outputs and evidence

Integrated runs write these files under `<outdir>/classification/mmseqs/`:

- `mmseqs_rescue.tsv`: query-to-reference alignments.
- `mmseqs_rescue.reverse.tsv`: reference-to-query alignments.
- `mmseqs_rescue.log`: the actual commands and full MMseqs output for each direction.
- `mmseqs_rescue.checkpoint.json`: completion state, input signatures and output hashes.
- `mmseqs_rescue.report.json`: query/reference paths and SHA256 hashes, parameters,
  algorithm version, cache status and search/evidence counts.

The main classification log distinguishes forward and reverse filtered
alignments, reciprocal pairs and families receiving rescue evidence. It also
reports the number of final class/order/superfamily decisions changed by rescue,
computed by rerunning the decision rules with rescue fields cleared while
holding all other evidence fixed. `previously_unknown_rescued` counts the subset
that changed from Unknown to a classified TE. These are decision-impact counts,
not biological accuracy measures or counts of genomic TE copies.

Per-family evidence remains in `<species>.evidence.tsv`, including
`rescue_target`, `rescue_order`, `rescue_superfamily` and alignment scores. The
final labels are in `<species>.classifications.tsv` and the classified FASTA.

## Cache validation and older assisted runs

Both input FASTA contents, search/filter parameters, executable and rescue
algorithm version are part of the signed search cache. Both directions must
complete and their output hashes must match before reuse. A valid empty search
can be reused. Changing CPU count alone does not change the evidence signature.
An isolated existing `mmseqs_rescue.tsv` from the older one-direction algorithm
is not accepted as a completed reciprocal cache; both searches are regenerated.
Failed or interrupted searches cannot authorize complete results.

Update the checkout, then rerun the original assisted command in the same output
directory. The integrated pipeline refreshes affected classification outputs
and downstream annotation products while preserving completed discovery,
extension and RepeatClassifier checkpoints. Existing compatible domain and Dfam
evidence is reused. A genuine unrelated biological-input/configuration change
still needs its normal checkpoint validation; the rescue upgrade does not
silently approve it. Allow the run to finish before comparing the resulting
annotations and inspecting the decision-impact counters.

Refinement now records the content hash of its raw RepeatMasker input. A new
annotation therefore invalidates old refinement results even when refinement
parameters and coordinates are unchanged. Older refinement manifests without
that hash are regenerated once; summaries and family inspection then refresh
from the rewritten refinement products. No manual directory deletion is needed
for this rescue upgrade.

## Standalone classifier

The same reciprocal path is available with a consensus FASTA:

```bash
drayte-classify --fasta TE_library.fa \
  --mmseqs-rescue-db Heliconius.classification.fa \
  --mmseqs-outdir classification_mmseqs \
  --output classifications.tsv --evidence-output evidence.tsv
```

`--mmseqs-rescue-tsv` remains the legacy interface for supplied within-library
evidence. It cannot be combined with `--mmseqs-rescue-db`: an unverified TSV
cannot replace the reference search's signed two-direction cache. A reference
search requires `--fasta`; the `--input` evidence-table mode has no query
sequences to search.
