# Preflight checks

`drayte` automatically checks the configured workflow after loading the YAML and
applying command-line overrides. These checks precede output-directory creation,
the output lock, archive restoration, the run log, and computational stages.
An enabled HELIANO stage with a missing executable at `structure.heliano.bin`
therefore fails before discovery or consensus extension begins.

To inspect the setup and then exit, run:

```bash
drayte --config config.yaml --check-only
```

Use the same overrides as the intended run. For example:

```bash
drayte --config config.yaml \
  --repeatmodeler-library /path/to/species-families.fa \
  --mmseqs-rescue-db /path/to/reference.classification.fa \
  --check-only
```

With `--check-only`, a successful check prints the selected executable and input
paths. Errors list the missing or inaccessible tools and files together, with
the corresponding configuration keys. Fix those settings or install the required
software in the environment used to run `drayte`, then repeat the check. Warnings also appear in
the report. The command returns a nonzero exit status when required checks fail.

The check reads files and locates commands without launching external tools.
It leaves the run directory and stored archives untouched. It creates no FASTA
indexes, pressed HMM indexes, DIAMOND databases, or FamDB exports.

## Requirements selected by the configuration

| Configuration or workflow | Checks |
| --- | --- |
| `genome`, `outdir` | A readable, nonempty genome file, an absolute output path, and access to the output directory or its existing writable parent. |
| Core pipeline | `perl`, `makeblastdb`, `blastn`, and `faToTwoBit` on `PATH`, plus `repeatmasker_bin`, `repeatclassifier_bin`, and `diamond_bin`. |
| Discovery without a supplied library | `BuildDatabase` and `RepeatModeler` under `repeatmodeler_dir`, plus `RepeatScout` and `build_lmer_table` under `repeatscout_dir`. |
| `repeatmodeler_library` or `--repeatmodeler-library` | The supplied library must be readable and nonempty. Discovery executables are skipped. Extension still requires the RepeatModeler utilities and scoring matrices below, and the `repeatscout_dir` setting remains required by the discovery wrapper. |
| Consensus extension | `repeatmodeler_extend_script`, the utilities and matrices under `repeatmodeler_dir`, `RAMExtend` under `repeat_after_me_dir`, and `cross_match` under `phrap_dir`. Directory and script paths used after changing the working directory must be absolute. |
| `ucsc_tools_dir` | Availability of the directory's UCSC helpers is advisory. Extension also requires `faToTwoBit` on `PATH`. RepeatModeler's `Linup` reads its helper location from its own `RepModelConfig`. |
| `repeatpeps_db_dir` | A readable, nonempty `RepeatPeps.lib`, including when a cached DIAMOND `.dmnd` database is present. The curation stage builds the DIAMOND database if it is absent. |
| `structure.heliano.enabled: true` | `structure.heliano.bin`, its `mmseqs_bin`, and HELIANO's commands on `PATH`, including `bedtools`, `rnabob`, `cd-hit-est`, `mafft`, `hmmsearch`, `getorf`, `gt`, `dialign2-2`, `Rscript`, and `bash`. |
| `classification.pfam_db` | The HMM file, all four pressed indexes (`.h3f`, `.h3i`, `.h3m`, `.h3p`), and `classification.hmmscan_bin`. Prepare indexes with `hmmpress` before checking. |
| Enabled `classification.dfam` | `nhmmer_bin` and the selected database. A FamDB directory also requires `famdb_bin`, a taxon, and readable, nonempty, uncompressed `*.h5` partitions. An exported HMM needs `taxon: null`. |
| `classification.mmseqs_rescue_db` or `--mmseqs-rescue-db` | A readable, nonempty reference FASTA and `classification.mmseqs_bin`. |

Disabled HELIANO adds no HELIANO requirements. Core ORF handling uses the internal
caller, so `getorf` and `cd-hit-est` enter the checks only through enabled
HELIANO. A disabled nested `classification.dfam` block adds no Dfam requirements
and supersedes the legacy flat `classification.dfam_db` setting. With no nested
block, a configured legacy `dfam_db` requires its HMM and `nhmmer`.

Unavailable Dfam commands and missing reference paths follow
`classification.dfam.on_missing: warn` when configured. Invalid Dfam settings,
such as supplying a single HDF5 partition instead of its directory, remain
errors. The historical `classification.enabled`, `preliminary`, and `final`
settings do not disable the integrated classification stage. Its configured
references determine the evidence checks.

YAML path strings do not expand shell variables. Replace placeholders such as
`CONDAENV` with actual paths in the run YAML. Executables used by stages that
change directory must use absolute paths or resolve through an absolute `PATH`
entry. Use absolute directory and script paths for extension. The input-library
and MMseqs rescue-reference settings support `~` expansion through their normal
runtime path handling.

## Resume and validation limits

Resumed runs receive the same configuration checks as new runs. Existing cache
files do not establish stage completion or remove tool requirements from
preflight. After the check succeeds, the normal pipeline restores managed
archives and validates each stage's checkpoints, inputs, and outputs.

Preflight verifies file accessibility and executable selection. Third-party
private configurations, shared libraries, software versions, reference contents
and taxonomic completeness, filesystem quota, and required R packages still need
their own checks. A readable HDF5 partition does not establish a complete Dfam
installation. An executable command can still fail when launched.

Python packages needed to import the runner must already be installed. Missing
packages can raise an import error before the preflight report is available.
Integrated summary generation uses pandas, matplotlib, and `tabulate` through
pandas' Markdown output. HELIANO has its own R dependencies. Use the
[environment verification instructions](environment.md#verify-the-installed-software)
and [Dfam setup instructions](environment.md#install-or-reuse-dfam-40) alongside
this check.
