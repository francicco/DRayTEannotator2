# Cluster environment and Dfam 4.0 references

The environment targets Linux/HPC. The release audit verified upstream
RepeatMasker 4.2.4, RepeatModeler 2.0.9 and FamDB 3.0.0. These versions are pinned
in `environment.yml`; HMMER is pinned to 3.4. RepeatModeler 2.0.9 explicitly
supports RepeatMasker 4.2.4 and Dfam 4.0.

Sources:
[RepeatMasker releases](https://www.repeatmasker.org/RepeatMasker/),
[RepeatModeler releases](https://github.com/Dfam-consortium/RepeatModeler/releases),
[RepeatMasker recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/repeatmasker/meta.yaml),
[RepeatModeler recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/repeatmodeler/meta.yaml),
[FamDB recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/famdb/meta.yaml).

## Create the environment

From the repository checkout:

```bash
micromamba create -n drayann-dfam40 -f environment.yml --strict-channel-priority
micromamba activate drayann-dfam40
python -I -m pip check
```

This creates a separate environment and leaves `drayann731` available for your
existing runs. The file uses conda-forge/bioconda and excludes implicit defaults.
The editable pip installation must run from this checkout and includes the FamDB
Python extra. A dry run, if desired, uses the same command with `--dry-run`.

Verify the software actually selected by this environment:

```bash
command -v python RepeatMasker RepeatModeler RepeatClassifier famdb.py nhmmer
RepeatMasker -h 2>&1 | sed -n '1,8p'
RepeatModeler -version
famdb.py --help
nhmmer -h
python -c 'import Bio, yaml, pandas, numpy, matplotlib, pyfaidx, pybedtools, h5py, tabulate; print("Python imports OK")'
perl -MGD::Simple -e 'print "GD::Simple OK\n"'
Rscript -e 'library(tidyverse); library(ggplot2); library(scales); library(Biostrings); library(BSgenome); cat("R libraries OK\n")'
```

The version command is tool-specific: RepeatModeler supports `-version`, whereas
RepeatMasker may reject that option and print usage. Use its help banner and the
installed Conda metadata rather than treating a rejected flag as an installation
failure.

If `pip check` reports unrelated Jupyter, scVI or other packages, first use the
exact active environment's Python in isolated mode:

```bash
"$CONDA_PREFIX/bin/python" -I -m pip check
```

`-I` excludes the user site and Python environment-variable paths; it does not
remove or install packages. If isolation clears the errors, they came from
outside this environment. If they persist, inspect the packages installed inside
the environment before changing dependencies. Do not add unrelated notebook
packages or downgrade NumPy solely to repair external package warnings.

For wrong tool banners or a `famdb.py` found outside the environment, inspect
both command selection and installed package metadata:

```bash
type -a python RepeatMasker RepeatModeler famdb.py
"$CONDA_PREFIX/bin/python" -I - <<'PYVERSIONS'
import json
import sys
from pathlib import Path
wanted = {"repeatmasker", "repeatmodeler", "famdb", "hmmer"}
found = {}
for path in (Path(sys.prefix) / "conda-meta").glob("*.json"):
    record = json.loads(path.read_text())
    if record["name"] in wanted:
        found[record["name"]] = record["version"]
print("Python:", sys.executable)
for name in sorted(wanted):
    print(name, found.get(name, "NOT INSTALLED"))
PYVERSIONS
```

Old versions in the metadata mean the environment itself needs correction.
Expected package versions with older executable banners instead call for checking
aliases, executable symlinks, PATH and PERL5LIB. Also verify the checkout and
`environment.yml` used to create the environment; an activation prompt alone
does not establish which manifest was installed.

FamDB 3.0 may print help and return exit status 1; that behavior alone does not
indicate a failed installation. The `famdb.py info` checks below are the database
checks that must succeed.

Perl is pinned to 5.32.1 to match the current repeat-tool builds. The previous
broad `perl-bioperl` metapackage was removed after the solver exposed its obsolete
Bio::Samtools/SAMtools dependency chain. Neither the current pipeline nor the
pinned upstream RepeatModeler source imports Bio::Perl modules; their local Perl
modules and declared package dependencies are retained.

The YAML now explicitly includes the pipeline's Python dependencies, HDF5 support,
HMMER, BLAST+, BEDTools, sequence tools, UCSC two-bit tools, RepeatAfterMe, MAFFT,
LTR tools and R libraries. RepeatModeler's remaining Perl dependencies are supplied
by its Bioconda package. `perl-gd` supplies the upstream PNG renderer's GD::Simple.
The extension script retains a custom Python renderer when installed, falling
back to upstream `visualizeAlignPNG.pl` when the custom file is absent.

## External tools and data are separate

This environment does not install the following reference data or external tools:

- Dfam 4.0 HDF5 partitions: see below.
- Pfam nucleotide-independent protein HMM profiles: configure `classification.pfam_db`
  and run `hmmpress` for the protein `hmmscan` database.
- `cross_match`: the customized consensus-extension script invokes this executable
  from `phrap_dir`. Keep your existing installation or install it separately.
- HELIANO, if `structure.heliano.enabled` is true. Keep the configured installation
  and its own dependencies. It is disabled in the default template.
- TE-Aid, if running that separate/manual inspection workflow. The current main
  curation runner copies existing extension artifacts; it does not invoke TE-Aid.

A package solve is not a full workflow test. Check these configured tools and run
the Hmel classification test before claiming cluster execution is validated.

## Install or reuse Dfam 4.0

`famdb=3.0.0` installs the query tool, not Dfam 4.0 database content. The official
[Dfam 4.0 README](https://dfam.org/releases/Dfam_4.0/families/FamDB/README.txt)
describes its version-3 format and compatibility with FamDB 3.0.0 / RepeatMasker
4.2.4.

For the curated-only Hmel test, the complete curated set listed by the official
[release directory](https://dfam.org/releases/Dfam_4.0/families/FamDB/) consists of:

```text
dfam40.0.h5
dfam40.curated.consensus.0.h5
dfam40.curated.hmm.0.h5
dfam40.curated.hmm.1.h5
```

The root holds release/taxonomy metadata and RepeatPeps; consensus models support
RMBlast/RepeatClassifier, while HMM models support Dfam/NHMMER validation. The
compressed HMM partition 1 alone is approximately 2 GB. This curated set is not
the complete curated-plus-uncurated Dfam release.

Reuse an already complete Dfam 4.0 directory where possible. To download the
curated set into a new directory, set a real writable path and run:

```bash
DFAM_DIR=/absolute/path/to/Dfam_4.0/FamDB
mkdir -p "$DFAM_DIR"
cd "$DFAM_DIR"
DFAM_URL=https://dfam.org/releases/Dfam_4.0/families/FamDB
for part in dfam40.0.h5 dfam40.curated.consensus.0.h5 dfam40.curated.hmm.0.h5 dfam40.curated.hmm.1.h5; do
  wget -c "$DFAM_URL/$part.gz" "$DFAM_URL/$part.gz.md5" || break
  md5sum -c "$part.gz.md5" || break
  if [ ! -f "$part" ]; then
    gzip -dk "$part.gz" || break
  fi
done
famdb.py -i "$DFAM_DIR" info
```

Confirm `Version : 4.0` and that all required components are present before
continuing. Additional uncurated files are needed if you enable uncurated models;
select those from the installed taxonomy/partition information rather than
assuming a partition number from older Dfam releases.

## Make RepeatClassifier use the same reference

The new DRayTE validation stage passes its `database` directory explicitly.
RepeatClassifier 2.0.9, however, calls its configured FamDB CLI without `-i`.
Therefore setting only `classification.dfam.database` does not configure the
reference used by RepeatClassifier.

For the Conda FamDB 3.0 installation, set its documented `famdb.conf` default data
directory to the same Dfam 4.0 directory. After setting `DFAM_DIR` above:

```bash
python - "$DFAM_DIR" <<'PYCONFIG'
import configparser
from pathlib import Path
import shutil
import sys
famdb_script = shutil.which('famdb.py')
if not famdb_script:
    raise SystemExit('famdb.py is not on PATH')
source = Path(sys.argv[1]).resolve()
if not (source / 'dfam40.0.h5').is_file():
    raise SystemExit('Missing Dfam 4.0 root partition')
# Bioconda's bin/famdb.py is a shell wrapper, not a symlink.
# The config belongs beside the real Python source in share/famdb-3.0.0.
famdb_home = Path(sys.prefix) / 'share' / 'famdb-3.0.0'
if not (famdb_home / 'famdb_globals.py').is_file():
    famdb_home = Path(famdb_script).resolve().parent
if not (famdb_home / 'famdb_globals.py').is_file():
    raise SystemExit('Locate the real FamDB source directory before configuring famdb.conf')
config_path = famdb_home / 'famdb.conf'
cfg = configparser.ConfigParser()
if config_path.exists():
    cfg.read(config_path)
    backup = config_path.with_name('famdb.conf.before-dfam40')
    if not backup.exists():
        backup.write_bytes(config_path.read_bytes())
if not cfg.has_section('famdb'):
    cfg.add_section('famdb')
cfg.set('famdb', 'FAMDB_DATA_DIR', str(source))
with config_path.open('w') as handle:
    cfg.write(handle)
print(config_path)
PYCONFIG
famdb.py info
```

This affects the selected environment's FamDB default, not a separate older
software installation. Check the installed RepeatModeler/RepeatMasker
configuration if they were built outside this Conda environment. A Conda package
update can replace `famdb.conf`; recheck `famdb.py info` afterward.

DRayTE curation separately requires a `RepeatPeps.lib` FASTA. Export it from the
same root partition into the directory used by `repeatpeps_db_dir`:

```bash
REPEATPEPS_DIR=/absolute/path/to/Dfam_4.0/repeatpeps
mkdir -p "$REPEATPEPS_DIR"
famdb.py -i "$DFAM_DIR" repeat_peps > "$REPEATPEPS_DIR/RepeatPeps.lib"
test -s "$REPEATPEPS_DIR/RepeatPeps.lib"
```

## Avoid accidentally using the old binaries

Installing new packages does not override absolute paths in a run YAML. In
particular, the previous default pointed at RepeatModeler 2.0.4. The template now
uses an explicit placeholder for the 2.0.9 directory and the `RepeatMasker`
command on PATH.

Inspect the active environment's paths:

```bash
python - <<'PYPATHS'
from pathlib import Path
import shutil
for name in ('RepeatModeler', 'RepeatMasker', 'RepeatClassifier', 'RepeatScout', 'RAMExtend', 'faToTwoBit', 'famdb.py', 'nhmmer', 'mmseqs'):
    path = shutil.which(name)
    print(f'{name}: {Path(path).resolve() if path else "MISSING"}')
PYPATHS
```

Set `repeatmodeler_dir` to the resolved RepeatModeler directory (Bioconda uses
`$CONDA_PREFIX/share/RepeatModeler`), `repeatmasker_bin` to the active RepeatMasker,
and `repeatclassifier_bin` to the active RepeatClassifier. For the packaged
RepeatScout, UCSC and RepeatAfterMe programs, their directory is normally
`$CONDA_PREFIX/bin`; put the actual expanded paths in the YAML. The config loader
does not expand shell variables in strings.

`repeatmodeler_extend_script` must point to this repository's customized
`davidExtendConsRAM.pl`, and `phrap_dir` to the directory containing your existing
`cross_match`. Set the new environment's `mmseqs` and HMMER paths where the run
configuration currently uses older absolute paths.

For Hmel, set `classification.dfam.database` to the verified directory and
`classification.dfam.taxon: Lepidoptera`. See the
[standalone/integrated test procedure](dfam-validation.md#first-cluster-test-on-an-existing-te-library).
Use a copy of your Hmel run YAML; changing tool/reference versions can invalidate
existing checkpoints, so do not bypass their input-mismatch checks.

## Local audit status

The environment declaration was checked against pipeline imports, invoked tools,
upstream source and current Bioconda recipes. After resolving an index-fetch
failure and removing the incompatible broad BioPerl metapackage, a strict-priority
Linux x86-64 Conda dry run resolved 504 packages. Cross-platform resolution used
Linux 4.18 / glibc 2.28 virtual-system settings, not a measured BluePebble node.
The solution included RepeatMasker 4.2.4, RepeatModeler 2.0.9, FamDB 3.0.0,
HMMER 3.4, Perl 5.32.1, GD 2.84 and h5py 3.16.0.

Run `uname -m` and `getconf GNU_LIBC_VERSION` on the target node, and perform the
Micromamba dry run there before installation. This local dry run validated Conda
package resolution; it did not install the environment, execute the editable pip
step, install databases or run the cluster genome workflow. Repository checks
passed with 487 tests passed / 3 skipped, and `perl -c davidExtendConsRAM.pl`
passed. Renderer dispatch tests cover custom Python and upstream Perl selection,
not full alignment rendering with an installed RepeatModeler stack.
