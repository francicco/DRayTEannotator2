# Recover a completed RepeatClassifier result from curation

Use this procedure when an old successful run reached curation, but later
retries rewrote the merged input timestamp or replaced the classified FASTA
before a modern completion checkpoint was created. A `running` marker records
an attempt, not proof that the process is currently alive or successfully done.

The usual resume path deliberately refuses to infer successful classification
from the presence of a FASTA alone. RepeatClassifier can orient consensuses by
reverse complementation; family identifiers and lengths alone cannot prove that
the classification process finished.

## Recovery source

This pipeline's curation stage writes **all** input records to
`curation/Final.RepeatModeler.Lib.fa`. It normalizes some classification headers
and strips descriptions, but does not discard or reorient the sequences. The
saved `curation/prioritize/family_table.tsv` contains each original classified
header and sequence length. Together with the completion manifests and unchanged
extension sources, these can support a reconstruction of the completed result.

The rewritten merged-input timestamp is not used as historical proof in this
procedure. Instead, the raw extension sources must predate the old successful
receipt, and their concatenated content must match the current merged input.

This is a recorded reconstruction from downstream evidence. It is not a claim
that RepeatClassifier ran again or that the restored FASTA is byte-identical to
the historical file. Descriptions and FASTA serialization may differ.

## Procedure

Use the original configuration for the same output directory. Preview while
inspecting the run; apply only after the current Mmex pipeline/worker attempt has
stopped. Other species can continue running in their own directories.

```bash
# Read-only: inspect all evidence and report eligibility/errors.
python -m drayte.pipeline.reclassify_recovery --config Mmex.drayannConfig.yaml

# Apply to an idle run, only when the preview is eligible.
python -m drayte.pipeline.reclassify_recovery --config Mmex.drayannConfig.yaml --apply

# Resume normally; the rebuilt completion checkpoint is checked as usual.
drayte --config Mmex.drayannConfig.yaml
```

The apply command acquires the same output lock as the pipeline. It cannot
exclude old jobs or external tools that never acquired that lock, so confirm
those workers are stopped as well. Do not start another concurrent retry or
manually copy files over the classified result.

## Required evidence

- A successful reclassification manifest identifying the expected merged and
  classified paths.
- A successful curation manifest identifying the curated library and family
  table in this same output tree.
- The current merged input exactly matches the selected extension-consensus
  files. Their modification/change times support their presence before the
  historical reclassification completed.
- The curated library and table have complete, unique family identifiers and
  unchanged lengths. Each curated sequence matches the current merged sequence
  either directly or by reverse complementation; unrelated substitutions are
  rejected.
- Original headers in the family table reproduce the normalized headers of the
  curated library. All source artifacts remain stable during validation.
- A mismatched modern **completed** checkpoint is not overridden. If a valid
  completed checkpoint already exists, recovery is unnecessary.

Missing, inconsistent, recently modified, or ambiguous evidence causes refusal
before replacement. The command does not force completion by editing timestamps
or setting a marker around an arbitrary current FASTA.

## Preserved evidence and provenance

Existing result/checkpoint/manifest files are backed up under
`reclassify/recovery_backups/` before replacement. The curated source and family
table remain unchanged. The recovery provenance record names the source files,
their hashes, the restored record counts, and the backup location. A normal
content-hashed completion checkpoint is written only after the reconstructed
FASTA validates.

This legacy procedure relies on saved workflow receipts and filesystem
chronology because historical content hashes were not recorded. It cannot
retroactively establish the exact old executable/database versions. New normal
runs record content hashes and require no such reconstruction.

Controlled interruptions of new RepeatClassifier attempts now trigger cleanup
of their isolated worker process group while leaving an incomplete checkpoint.
This change cannot stop workers already launched by older code.
