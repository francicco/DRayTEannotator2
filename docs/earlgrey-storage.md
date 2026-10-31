# EarlGrey storage

This repository provides a separate, reversible storage command for EarlGrey
trees. Native EarlGrey does not know about these archives. This command never
modifies EarlGrey's installed source or changes TE detection/classification.
The layout and completion checks follow the
[official EarlGrey pipeline](https://github.com/TobyBaril/EarlGrey/blob/main/earlGrey).
Older or different layouts that do not satisfy these checks are refused.

## Existing runs

Use the Python environment where DrayTE is installed. Preview one outer output
directory (for example `Avcr_Earlgrey`) or its inner `Avcr_EarlGrey` directory:

```bash
python -m drayte.storage.earlgrey --root Avcr_Earlgrey
```

Preview changes no output files. It reports `eligible`, validation errors, an
archive plan and current allocated bytes/file counts. Completion requires the
native final-stage log receipt, final annotation, library and divergence/summary
products, an independent discovery library, and the raw RepeatMasker hit/table
pair. Exported library and BED/GFF contents must match their stage sources.
Existing archives are verified and symlink dependencies are checked.

An 11–17 KB output tree can be a precreated or unfinished set of stage
directories. Small size does not establish a successful run, and final summaries
alone can occupy hundreds of megabytes. There is no universal target size.

From the directory holding the species outputs, preview them all:

```bash
for directory in *_Earlgrey; do
  [ -d "$directory" ] || continue
  python -m drayte.storage.earlgrey --root "$directory"
done
```

For an eligible run, first confirm that no job or worker is using that specific
output tree. The `--idle` flag records this assertion; a new lock cannot detect
older native EarlGrey processes. Then apply:

```bash
python -m drayte.storage.earlgrey --root Avcr_Earlgrey --apply --idle
```

Do not apply indiscriminately to every species while jobs are running. Other
species can continue running if they do not read/write the selected tree.

## Retained and archived data

The following stay directly available:

- The complete `*_summaryFiles` directory, including BED/GFF, library and plots.
- The complete `*_Database` directory, including `*-families.fa` used as a DrayTE input.
- Top-level strained libraries, native logs, and completion/skip markers.
- The final raw RepeatMasker `.out` and `.tbl` used in benchmarking, and the
  already compressed `.cat.gz`.
- HELIANO's top-level `RC.representative.*` outputs.

Only known bulky paths are archived: RepeatModeler `RM_*` scratch, strainer
`TS_*` work, merge intermediates whose final BED/GFF exports match, landscape
work, final RepeatMasker `.align`/`.masked`, and HLE1/HLE2 or MMseqs
work beneath completed HELIANO runs. Unknown directories stay untouched.
Retained links into a candidate block compaction. External symlinks in an
archive source also block it. No new EarlGrey policy permanently discards data.

The command writes `.storage/earlgrey_report.json` in the inner run directory,
with actual before/after allocated bytes and file counts. Verified compressed
payloads and the recovery catalog live alongside it under `.storage/`.
Never remove this directory: it is part of the resumable dataset.

## Future jobs and resume

Use this wrapper around the original EarlGrey command in the job script:

```bash
python -m drayte.storage.earlgrey --root Avcr_Earlgrey --idle --run \
  earlGrey -g genome.fa -s Avcr -o Avcr_Earlgrey -t 16
```

Retain all original EarlGrey options. The species and outer `-o` directory must
match the wrapper's selected root. The wrapper acquires its lock, verifies and
restores managed files to their original paths, runs EarlGrey, then compacts
only if EarlGrey exits successfully and completed outputs validate. An existing
nonempty run requires `--idle`; a new run does not. A failed native command
leaves restored working files available for investigation and retry.

Before using native EarlGrey directly, restore explicitly:

```bash
python -m drayte.storage.earlgrey --root Avcr_Earlgrey --restore --apply --idle
```

Restoration verifies checksums and preserves original paths and modification
times. Modified loose files are never overwritten. Corrupt/missing archives
stop restoration. Archive creation needs temporary quota headroom while
originals remain; failure before verification leaves those originals intact.
Resume expands archived files, so working space can exceed compact idle space.
Already compressed files may save few bytes; gains depend on actual data.

The module form works after a repository update without reinstalling. Installing
the updated package also registers the equivalent `earlgrey-storage` shortcut.
