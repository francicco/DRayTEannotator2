# Storage management

This page describes storage management integrated into DrayAnn. EarlGrey output
has a separate [single-file utility](earlgrey-storage.md) that runs independently
of DrayAnn.

Storage has two different costs: bytes and filesystem objects. One large
alignment consumes disk space; thousands of per-family work files also consume
inode quota and make backups and directory scans expensive. Compact mode reduces
both while preserving evidence needed to resume.

The final delivery folder is not a measure of a pipeline's total footprint.
For example, the supplied Hmel directory listing showed 2.7 GiB across EarlGrey's full
run tree even though a particular final folder can be small. Compare allocated
bytes and file counts across the whole run, including retained work directories.

## New runs

Add to the existing pipeline YAML:

```yaml
storage:
  mode: compact
```

Or use the command-line override:

```bash
drayte --config Eisa.config.yaml --storage-mode compact
```

`keep` is the compatibility default. It retains loose products. Both modes
automatically restore any previously managed archives before running stages.
Compaction is opt-in; no existing cluster result is changed just by updating
the repository.

## Lifecycle

1. Acquire an exclusive lock for the output directory. An updated pipeline and
   the storage command cannot modify the same run concurrently.
2. Restore and verify any managed archives **before** discovery, extension,
   curation, or classification checks. Original paths and modification times
   are restored, so archived evidence is not mistaken for missing input.
3. Run/resume the normal stages. In compact mode, after successful discovery,
   delete redundant `discovery/rmodeler_dir/RM_*` work directories only when an
   independent top-level library and RepeatMasker hit/table outputs exist.
   Configured inputs and retained symlinks pointing into a candidate prevent
   deletion. After completed HELIANO/classification, their MMseqs scratch can
   likewise be removed. These early removals reduce peak storage during a run.
   Once extension completes with no failed families, archive its entire work
   tree while retaining regular files needed by postprocessing, curation,
   classification and inspection. This avoids carrying all extension scratch
   through HELIANO, classification and final RepeatMasker. It does not reduce
   the space needed while extension itself is running.
4. After the **whole pipeline succeeds**, validate retained products, extension
   summaries, RepeatClassifier completion, classification biological hashes,
   and RepeatMasker completion records, then archive the remaining bulky evidence.
   An unfinished stage's evidence is kept loose; a completed extension stage
   can be compacted while later stages are still unfinished.
5. Create compressed archives, read them back and verify every file checksum,
   then commit the recovery catalog **before** removing the loose originals.
   Never remove originals merely because a compression command returned zero.
6. Write `.storage/storage_report.json` with actions, archive locations, and
   before/after allocated-byte and file counts. Early scratch cleanup is logged
   in `.storage/scratch_actions.jsonl` and `drayte.log`.
   Early extension savings are measured separately in `.storage/extension_report.json`.

For existing runs compacted through the separate DrayTE storage command,
original CLI-only inputs may be absent from the YAML. Scratch is therefore archived reversibly
rather than deleted. A resumed pipeline restores it before checking the original
`--repeatmodeler-library`/reference paths. Automatic scratch deletion is limited
to pipeline runs where the effective runtime inputs are known.

Failed computational stages retain their unfinished evidence for retry. An
interruption during archive removal leaves a verified archive and recovery
catalog. Restoration fills in missing files, refuses conflicting modified
files, and never silently overwrites them. Corrupt or missing archives stop
resume explicitly.

## What stays loose

- The original genome input and normalized working genome.
- Discovery and final TE libraries, including the classified FASTA.
- Discovery and final RepeatMasker `.out/.tbl`, and final `.cat.gz` (already
  compressed). Existing GFF output remains loose.
- Refined and filtered BED/GFF/TSV annotations and summaries.
- Stage manifests, completed extension summary, final extension consensuses,
  and the cached extension-query FASTA under `extension/blastfiles`.
- During execution, extension work retains normalized consensus/alignment/image
  files, raw postprocessing sources (including `repseq.unextended`, which carries
  hit counts), and `repam-ranges.tsv`/`repam-repseq.fa`. Configured file inputs and
  files linked from elsewhere in the output tree are also retained. The full
  original tree is recorded in its verified archive. A referenced directory or
  symlinked required file blocks this early compaction instead of breaking it.
- Merged Pfam/Dfam evidence and family classification/evidence tables.
- Existing portable inspection `.tgz` deliverables from older runs. New compact
  runs avoid creating these duplicate dereferenced copies; managed archives
  already preserve the family evidence.

## What is archived

- `extension/extensionwork`, `extract_align`, `extendlogs`, and
  `images_and_alignments`.
- `curation/te-aid`.
- Pfam `chunks`/`hmmscan_chunks` and Dfam `chunks`/`tblouts`.
- Discovery masked FASTA, final masked FASTA, and final RepeatMasker alignment.
- Loose family-inspection category directories.

These are reversible archives, not discarded scientific evidence. The
family-inspection index still records original loose directory paths; use
`--restore --apply` before opening those directories or their symlinks. The
storage catalog/report maps each archived path to its payload. Symlinks within
the run are preserved rather than duplicating their target contents. External
symlinks in archive sources are rejected and require review.

## Existing completed runs

Use the same configuration and tool settings that produced the run. Update the
local checkout, then use the existing Python environment. The module command
works without recreating the executable entry point:

```bash
# Preview only: show candidates and validation result; modify no output files.
python -m drayte.storage.cli --config Eisa.config.yaml

# Apply only once you have confirmed no job is using this output directory.
python -m drayte.storage.cli --config Eisa.config.yaml --apply

# Restore all loose evidence for inspection or non-pipeline tools.
python -m drayte.storage.cli --config Eisa.config.yaml --restore --apply
```

To register the `drayte-storage` shortcut after updating an editable installation:

```bash
python -m pip install --no-deps --no-build-isolation -e .
```

The CLI deliberately processes one configuration/output directory at a time.
Preview each species before applying a batch. It only operates on DrayAnn2's
own output tree; external EarlGrey results and externally supplied libraries
are not compacted or removed by this command. For EarlGrey trees, copy the
standalone [`compact_earlgrey.py`](../compact_earlgrey.py) file anywhere and follow
the [EarlGrey storage instructions](earlgrey-storage.md). That utility requires
only standard-library Python 3.10+ and no DrayAnn installation.

Old running jobs do not hold the new lock. The lock cannot prove that a legacy
job, an external tool, or a person is not reading the directory. Confirm those
jobs are idle before using `--apply`. No live cluster directories were altered
during development.

## RepeatMasker resume records

New successful runs record hashes of the genome, library, biological command
options, and required output contents. A marker is written as `running` before
execution and becomes `complete` only after successful output validation.
An isolated `.tbl` from a quota failure cannot satisfy discovery completion.
Changing the timestamp of an otherwise identical library does not rerun a
completed modern annotation.

Legacy results have no historical input hashes. Complete products may be
adopted only when the required outputs exist and the last recorded attempt is
not failed or unfinished; final annotation also requires its library not to be
newer than its hits. Compaction labels this compatibility adoption explicitly
in the new checkpoint. It cannot retroactively establish historical tool or
database provenance. Ambiguous/incomplete legacy runs are refused rather than
being silently declared complete.

## Capacity and recovery

- Compression needs temporary free space for an archive while its originals
  still exist. During a compact pipeline run, completed disposable scratch is
  pruned first to make headroom. Configuration-based legacy compaction instead
  preserves scratch in archives, so it may require additional headroom. If archive
  creation fails, originals remain. You may still need another filesystem to
  bootstrap cleanup of a completely full quota.
- Resume temporarily expands archived evidence. This version reduces idle
  storage and prunes proven scratch during execution; it does not promise that
  peak working storage equals compacted storage.
- Do not move a run to a new absolute path as part of cleanup. Existing manifests
  and symlinks may retain absolute paths. Relocation is a separate operation.
- Do not delete `.storage/archives.json` or its archive payloads. They are part
  of the resumable dataset, not temporary files.
- `storage_report.json` is a snapshot of the last compaction. The live archive
  catalog is `archives.json`; after a successful restore its consumed payloads
  and records are removed to avoid retaining a second copy.
- If compaction is interrupted, use `--restore --apply` (or rerun `drayte`) before
  attempting another compaction. Never force-overwrite modified loose evidence.
- Actual savings depend on the genome and evidence. Tests verify round trips,
  recovery after interrupted removal, corrupt-archive refusal, preservation of
  checkpoint fingerprints, and RepeatMasker reuse after restore. No full-size
  Eisa/Hmel compaction has yet been measured on the cluster.
