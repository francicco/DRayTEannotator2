# Reference-assisted classification

Use a classified TE FASTA as homology evidence for the families being classified:

```bash
drayte --config Hmel.drayannConfig.Assisted.yaml \
  --repeatmodeler-library "$RModelerFasta" \
  --mmseqs-rescue-db "${HELICONIINILIB%.lib}.classification.fa"
```

The reference flag overrides `classification.mmseqs_rescue_db` in the configuration.
Reference headers need RepeatMasker labels, for example
`>Helitron-1_Hmel#RC/Helitron`. The reference supplies classification evidence;
its sequences are not appended to the final annotation library. Protein-domain,
Dfam, structural and header evidence still participates in the classification
rules. A changed label measures a classifier decision, not independently
established biological accuracy.

## One-way label transfer

External-reference rescue uses one-way matching: classification-query FASTA
against the curated reference. Each query is assessed independently, so several
predicted consensuses can receive classification evidence from the same
reference family. Query and reference identifiers remain separate, including when
the same clean family name occurs in both libraries.

The search prefilter remains at 75% identity and 50% coverage of both sequences. Retained alignment
filters remain at 80% identity, 500 aligned bases, 40% query coverage and 30%
reference coverage. These filters were kept unchanged to isolate the effect of
the transfer policy. A short or partial match can still fail these filters even
when it represents real homology. The integrated pipeline can override them
with `classification.mmseqs_min_identity`, `mmseqs_min_aln_len`,
`mmseqs_min_query_cov` and `mmseqs_min_target_cov`.
Lowering a retained-hit threshold cannot recover matches already excluded by
the search prefilter.

For each query, the policy examines the best qualified match to each reference
and all competing references scoring at least 95% of the best bit score. This
fraction is configurable with `--mmseqs-rescue-score-fraction`, or
`classification.mmseqs_rescue_score_fraction` in the integrated configuration.
Its valid range is greater than zero and at most one. The 95% window is a
configurable competition heuristic, not a validated biological family boundary.

- Competitive references agreeing on class and order support that assignment.
- Conflicting class/order labels, or a competitive reference without a usable
  classification, cause the policy to abstain.
- References agreeing on order but differing in superfamily support only the
  broader order assignment.
- Superfamily transfer additionally requires agreement among the competitive
  references and at least 80% identity and 80% coverage of both sequences.

Matching a curated family supplies evidence about a query's class, order and
superfamily. It does not establish that the query and reference are the same
biological TE family. A local match can cover a shared domain or only part of a
consensus. Exact family delimitation needs its own criteria and validation.

## Search outputs and evidence

Integrated runs write these files under `<outdir>/classification/mmseqs/`:

- `mmseqs_rescue.tsv`: query-to-reference search alignments, before retained-hit
  filtering.
- `mmseqs_rescue.log`: actual commands and full MMseqs output.
- `mmseqs_rescue.checkpoint.json`: completion state, input signatures and current
  output hashes.
- `mmseqs_rescue.report.json`: policy, query/reference paths and SHA256 hashes,
  parameters, algorithm version, cache status, search/evidence counts and
  per-query decisions.

The report's `query_decisions` object records each query's selected target,
competing labels, alignment metrics and reasons for rejected hits. A
`no_search_hit` reason means the search returned no alignment after its own
filters; it cannot establish an absence of biological homology.

Logs distinguish alignments, selected rescue evidence and final family labels
actually changed by that evidence. The change counts come from rerunning the
decision rules with rescue fields cleared while holding other evidence fixed.
`previously_unknown_rescued` counts the subset changing from Unknown to a
classified TE. These counts describe family decisions, not genomic TE copies
or biological accuracy.

Per-family evidence remains in `<species>.evidence.tsv`, including
`rescue_target`, `rescue_order`, `rescue_superfamily` and alignment scores. Final
labels are in `<species>.classifications.tsv` and the classified FASTA.

## Cache validation and older assisted runs

Input FASTA contents, search/filter settings, transfer policy, executable and
algorithm version are part of the signed search cache. Rescue requires
a completed, verified query-to-reference search.
Valid empty searches can be reused. Changing CPU count alone does not change
the evidence signature. An isolated TSV without a matching completed checkpoint
cannot authorize cached reference rescue. Failed or interrupted searches cannot
authorize complete results.

Completed version 2 one-way checkpoints remain reusable. Historical reciprocal
results are refreshed using the current one-way policy. Older reverse-search
files remain untouched and are never read or required. A legacy
`classification.mmseqs_rescue_mode: reciprocal` setting is ignored with a warning.
The previously supplied `--mmseqs-rescue-mode one-way` command is accepted for
compatibility; no selector is needed for current commands.

Update the checkout, then rerun against the same output directory:

```bash
git -C /user/work/tk19812/software/DRayTEannotator2 pull \
  --ff-only origin feature/parallel-consensus-extension

drayte --config Hmel.drayannConfig.Assisted.yaml \
  --repeatmodeler-library Hmel_Earlgrey/Hmel_EarlGrey/Hmel_Database/Hmel-families.fa \
  --mmseqs-rescue-db Heliconius.classification.fa
```

The pipeline refreshes affected classification and downstream annotation
products while preserving completed discovery, extension and RepeatClassifier
checkpoints. Compatible domain and Dfam evidence is reused. Unrelated changes
to biological inputs still receive normal checkpoint validation. Allow the run
to finish before comparing annotation coverage and decision-impact counts.

Refinement records the content hash of its raw RepeatMasker input. A new
annotation invalidates old refinement results even when refinement parameters
are unchanged. Older refinement manifests without that hash regenerate once;
summaries and family inspection then refresh from the rewritten products.
No manual directory deletion is needed for the rescue upgrade.

Using the same curated reference as both classifier input and evaluation
reference measures reference agreement. Use independent validation to estimate
accuracy, and report detection coverage separately from classification agreement.

## Standalone classifier

The same matching policy is available with a consensus FASTA:

```bash
drayte-classify --fasta TE_library.fa \
  --mmseqs-rescue-db Heliconius.classification.fa \
  --mmseqs-rescue-score-fraction 0.95 \
  --mmseqs-outdir classification_mmseqs \
  --output classifications.tsv --evidence-output evidence.tsv
```

`--mmseqs-rescue-tsv` retains the legacy within-library evidence interface. It
cannot be combined with `--mmseqs-rescue-db`: an unverified TSV cannot replace
the reference search's signed cache. A reference search requires `--fasta`;
the `--input` evidence-table mode has no query sequences to search.
