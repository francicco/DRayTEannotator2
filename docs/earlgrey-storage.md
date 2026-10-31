# EarlGrey storage

[`compact_earlgrey.py`](../compact_earlgrey.py) is a standalone, single-file
utility for reversibly reducing EarlGrey output trees. Copy this one file
anywhere and run it with Python 3.10 or newer. It uses only the Python standard
library: no DrayAnn installation, `pip` packages, or Conda environment is needed.
The optional `--run` wrapper still requires your existing EarlGrey installation
and its usual dependencies.

Native EarlGrey does not know about these archives. This utility never
modifies EarlGrey's installed source or changes TE detection/classification.
The layout and completion checks follow the
[official EarlGrey pipeline](https://github.com/TobyBaril/EarlGrey/blob/main/earlGrey).
Older or different layouts that do not satisfy these checks are refused.

## Existing runs

Replace `/path/to/compact_earlgrey.py` below with the location of your copied
script. Preview one outer output directory (for example `Avcr_Earlgrey`) or
its inner `Avcr_EarlGrey` directory:

```bash
python3 /path/to/compact_earlgrey.py --root Avcr_Earlgrey
```

If the root has a different name, add `--species Avcr` explicitly.

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
  python3 /path/to/compact_earlgrey.py --root "$directory"
done
```

For an eligible run, first confirm that no job or worker is using that specific
output tree. The `--idle` flag records this assertion; a new lock cannot detect
older native EarlGrey processes. Then apply:

```bash
python3 /path/to/compact_earlgrey.py --root Avcr_Earlgrey --apply --idle
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
python3 /path/to/compact_earlgrey.py --root Avcr_Earlgrey --idle --run \
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
python3 /path/to/compact_earlgrey.py --root Avcr_Earlgrey --restore --apply --idle
```

Restoration verifies checksums and preserves original paths and modification
times. Modified loose files are never overwritten. Corrupt/missing archives
stop restoration. Archive creation needs temporary quota headroom while
originals remain; failure before verification leaves those originals intact.
Resume expands archived files, so working space can exceed compact idle space.
Already compressed files may save few bytes; gains depend on actual data.

The standalone script uses the same archive format and recovery catalog as the
previous packaged utility. Existing managed archives can be previewed, verified,
and restored with this script; no conversion or DrayAnn installation is needed.
DrayAnn's integrated storage mode remains separate and unchanged.
