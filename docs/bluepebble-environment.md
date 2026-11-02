# Dedicated DrayAnn environment on BluePebble

`environment.bluepebble.yml` creates a native Linux environment named
`drayann731`. It has its own Python and external tools and does not require
activating the existing Python virtual environment under
`$APPS2/DRayTEannotator2/DRayTEannotator2/` or using the `earlgrey731`
environment at runtime.

The comparison toolchain pins RepeatMasker 4.2.4, RepeatModeler 2.0.9, FamDB
3.0.0, and Bioconda RMBlast 2.17.0 build 1. The latter includes the upstream
RMBlast 2.17.1 patch. Record both the package build and `rmblastn -version`
output when comparing environments. The specification also includes all
Python dependencies from `pyproject.toml`, RepeatAfterMe, HMMER, UCSC tools,
DIAMOND, MMseqs2, BEDTools, EMBOSS, and HELIANO. HELIANO brings its required R
runtime; the former broad tidyverse/Bioconductor dependency set is unnecessary
for the standard DrayAnn pipeline.

## Create and install

Run on BluePebble, from a shell that has Micromamba available. If the old Python
virtual environment is active, run `deactivate` before these commands. Keep the
existing environment and old analysis directories in place.

```bash
cd "$APPS2/DRayTEannotator2"
git pull --ff-only origin feature/parallel-consensus-extension

micromamba create --override-channels -c conda-forge -c bioconda \
  --strict-channel-priority -f environment.bluepebble.yml
micromamba activate drayann731

python -m pip install --no-deps --no-build-isolation -e .
python -m pip check
```

The editable installation uses this checkout, so later code updates should be
recorded with their Git commit. It does not install another copy of the TE
tools. Do not source the old `DRayTEannotator2/bin/activate` inside the new
environment.

## Verify the installation

```bash
command -v python drayte RepeatMasker RepeatModeler RepeatClassifier rmblastn
python -c 'import sys, h5py; print(sys.executable); print("h5py", h5py.__version__)'
RepeatMasker -version
RepeatModeler -version
rmblastn -version
drayte --help
drayte-classify --help
perl -MGD::Simple -e 'print "GD::Simple available\n"'

(
  for tool in BuildDatabase RepeatScout RAMExtend hmmscan nhmmer \
    mmseqs diamond cd-hit-est getorf bedtools blastn makeblastdb \
    faToTwoBit twoBitToFa twoBitInfo heliano; do
    command -v "$tool" || exit 1
  done
)

micromamba list -n drayann731 |
  grep -E '^[[:space:]]*(repeatmasker|repeatmodeler|rmblast|famdb)[[:space:]]'
micromamba list -n earlgrey731 |
  grep -E '^[[:space:]]*(repeatmasker|repeatmodeler|rmblast|famdb)[[:space:]]'
```

Executable paths for the new toolchain should be inside the new environment.
The editable `drayte` command should import this checkout using that
environment's Python. Compare the listed core package versions/builds with
`earlgrey731`; version pins alone do not prove that an older EarlGrey
installation selected the same RMBlast build.

Save installation provenance after a successful install:

```bash
mkdir -p environment-records/drayann731
micromamba list -n drayann731 --explicit --md5 \
  > environment-records/drayann731/conda-explicit.txt
python -m pip freeze > environment-records/drayann731/python-freeze.txt
git rev-parse HEAD > environment-records/drayann731/drayte-commit.txt
```

The explicit specification records exact Conda package versions, builds, and
checksums for Linux. The editable package and analysis inputs still need their
own recorded Git revision and file provenance.

## Configure paths for new comparison runs

An activated environment does not override absolute executable paths in YAML.
Existing configurations may still point to RepeatModeler 2.0.4 and a manually
installed RepeatMasker. Copy each analysis configuration to a new benchmark
configuration, give it a fresh output directory, and resolve these settings
against the activated `drayann731` environment:

| YAML setting | New environment location |
| --- | --- |
| `repeatmodeler_dir` | `share/RepeatModeler` under the environment prefix |
| `repeatmasker_bin` | `bin/RepeatMasker` |
| `repeatclassifier_bin` | `bin/RepeatClassifier` |
| `repeatscout_dir` | `bin` |
| `ucsc_tools_dir` | `bin` |
| `repeat_after_me_dir` | `bin` |
| `diamond_bin` | `bin/diamond` |
| `getorf_bin` | `bin/getorf` |
| `cdhit_est_bin` | `bin/cd-hit-est` |
| `structure.heliano.bin` | `bin/heliano` |
| `structure.heliano.mmseqs_bin` | `bin/mmseqs` |
| `classification.mmseqs_bin` | `bin/mmseqs` |
| `classification.hmmscan_bin` | `bin/hmmscan` |
| `classification.nhmmer_bin` | `bin/nhmmer` |

Print the activated prefix with `printf '%s\n' "$CONDA_PREFIX"` and use that
resolved absolute path in the copied YAML. YAML does not expand a literal
`$CONDA_PREFIX` variable.

Use the repository's custom `davidExtendConsRAM.pl` for
`repeatmodeler_extend_script`. RepeatModeler 2.0.9 supplies its own extension
utilities under `share/RepeatModeler/util`; the DrayAnn wrapper is a separate
file in this repository.

The custom extension step also requires the existing Phrap `cross_match`
executable. Keep `phrap_dir` pointed at its installed directory and verify it
is executable before launching. Phrap is not installed by this specification.

## Dfam and other databases are a separate setup step

Installing an environment does not install or migrate analysis databases.
Use the shared external Dfam 4.0 directory rather than downloading another
complete copy into each environment, and ensure both RepeatMasker and
RepeatModeler can find it through their FamDB configuration.

There are distinct database consumers:

- RepeatClassifier uses the FamDB root and curated consensus components. The
  root supplies RepeatPeps; the consensus component supplies nucleotide
  reference sequences.
- DrayAnn curation requires a plain `RepeatPeps.lib` file in the directory
  specified by `repeatpeps_db_dir`, with a DIAMOND index built there. That
  directory needs a Dfam 4.0 protein export when migrating the curation data.
- DrayAnn's extra HMM classification uses the plain HMM file given by
  `classification.dfam_db`. Root plus curated consensus files do not provide
  those profile HMMs. The HMM export, metadata mapping, and cache migration must
  be validated separately before describing this stage as Dfam 4.0 enabled.
- `classification.pfam_db` remains a separate protein-domain HMM database.

Only submit new comparison runs after these tool paths and database inputs
have been verified. A successful environment solve or `--help` invocation does
not establish that the complete biological pipeline has run successfully.

## FamDB cannot import h5py

`ModuleNotFoundError: No module named 'h5py'` from `famdb.py` means the Python
selected for FamDB lacks its HDF5 dependency. This can happen when the old
DrayAnn virtual environment is active while the YAML selects external tools
from `earlgrey731`. Changing `repeatmasker_bin` alone does not select the
Python that FamDB needs.

The dedicated environment includes `h5py`. Verify it without activating the
old virtual environment:

```bash
micromamba run -n drayann731 python -c \
  'import sys, h5py; print(sys.executable); print("h5py", h5py.__version__)'
```

Launch DrayAnn through `micromamba run -n drayann731 drayte` with the copied
configuration pointing to that same environment. If `drayann731` does not yet
exist, complete the creation and editable installation above first. A missing
Python module does not require downloading the Dfam data again.

The pipeline also now omits the previous bare `-gc` argument from both
RepeatMasker stages. RepeatMasker requires a value for `-gc`; without one it
can consume the following `-gff` argument. The corrected command uses
RepeatMasker's default matrix selection and retains `-gff`. It does not force
`-gccalc`. Checkpoints signed with the old argument list will be rerun.

## Sources

- [RepeatMasker Bioconda recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/repeatmasker/meta.yaml)
- [RepeatModeler Bioconda recipe and configured utility layout](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/repeatmodeler/build.sh)
- [RMBlast build and 2.17.1 patch](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/rmblast/meta.yaml)
- [Micromamba explicit specifications](https://mamba.readthedocs.io/en/stable/user_guide/micromamba.html#explicit-spec-files)
- [FamDB 3.0 documentation](https://github.com/Dfam-consortium/FamDB/blob/3.0.0/README.md)
- [RepeatMasker 4.2.4 option parser](https://github.com/Dfam-consortium/RepeatMasker/blob/3f49a049c7a7e843b2da1b708fe06bb39c90035b/RepeatMasker#L392)
