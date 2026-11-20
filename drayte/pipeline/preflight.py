"""Read-only checks of the tools and inputs selected by a pipeline configuration."""
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from .config import ConfigError, PipelineConfig


class PreflightError(ConfigError):
    """All setup problems found before the output tree is changed."""


@dataclass
class PreflightReport:
    tools: list[tuple[str, str, str]] = field(default_factory=list)
    inputs: list[tuple[str, str, str]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def format(self, *, detailed: bool = False) -> str:
        lines = [f"Pre-flight check passed ({len(self.tools)} tools, {len(self.inputs)} input paths)."]
        if detailed:
            for kind, checks in (("Tool", self.tools), ("Input", self.inputs)):
                lines.extend(f"  {kind}: {label}: {path} [{setting}]" for label, path, setting in checks)
        lines.extend(f"  WARNING: {warning}" for warning in self.warnings)
        return "\n".join(lines)


class _Checks:
    def __init__(self):
        self.report = PreflightReport()
        self.errors: list[str] = []

    def problem(self, message: str, *, warn: bool = False):
        (self.report.warnings if warn else self.errors).append(message)

    def text(self, value, setting: str, *, warn: bool = False) -> str | None:
        if not isinstance(value, str) or not value.strip():
            self.problem(f"{setting}: set a non-empty path or command in the YAML.", warn=warn)
            return None
        return value

    def executable(self, label: str, command, setting: str, *, warn: bool = False, changed_cwd: bool = False):
        command = self.text(command, setting, warn=warn)
        if command is None:
            return
        if changed_cwd and os.sep in command and not Path(command).is_absolute():
            self.problem(f"{label}: use an absolute executable path because this stage changes its working directory: {command} [{setting}].", warn=warn)
            return
        resolved = shutil.which(command)
        if resolved is None:
            detail = "not found on PATH" if os.sep not in command else "missing or not executable"
            self.problem(f"{label}: {detail}: {command} [{setting}].", warn=warn)
            return
        if changed_cwd and not Path(resolved).is_absolute():
            self.problem(f"{label}: resolves through a relative PATH entry: {resolved}; use an absolute PATH entry or executable path [{setting}].", warn=warn)
            return
        self.report.tools.append((label, str(Path(resolved).resolve()), setting))

    def file(self, label: str, value, setting: str, *, warn: bool = False, absolute: bool = False):
        value = self.text(value, setting, warn=warn)
        if value is None:
            return
        path = Path(value)
        if absolute and not path.is_absolute():
            self.problem(f"{label}: use an absolute path because extension changes its working directory: {path} [{setting}].", warn=warn)
            return
        try:
            if not path.is_file():
                raise ValueError("file is missing or is not a regular file")
            if path.stat().st_size == 0:
                raise ValueError("file is empty")
            with path.open("rb") as handle:
                handle.read(1)
        except (OSError, ValueError) as error:
            self.problem(f"{label}: {error}: {path} [{setting}].", warn=warn)
            return
        self.report.inputs.append((label, str(path.resolve()), setting))

    def directory(self, value, setting: str, *, absolute: bool = False) -> Path | None:
        value = self.text(value, setting)
        if value is None:
            return None
        path = Path(value)
        if absolute and not path.is_absolute():
            self.problem(f"Use an absolute directory path because extension changes its working directory: {path} [{setting}].")
            return None
        if not path.is_dir() or not os.access(path, os.R_OK | os.X_OK):
            self.problem(f"Directory is missing or unreadable: {path} [{setting}].")
            return None
        self.report.inputs.append(("directory", str(path.resolve()), setting))
        return path

    def mapping(self, value, setting: str) -> dict:
        if not isinstance(value, dict):
            self.problem(f"{setting}: expected a YAML mapping.")
            return {}
        return value

    def finish(self) -> PreflightReport:
        if self.errors:
            lines = [f"  - {error}" for error in self.errors]
            lines.extend(f"  WARNING: {warning}" for warning in self.report.warnings)
            lines.append("Fix the listed YAML settings or install the missing tools in the environment used to run drayte.")
            raise PreflightError("\n".join(lines))
        return self.report


def _check_dfam(checks: _Checks, classification: dict):
    if "dfam" not in classification:
        if classification.get("dfam_db"):
            checks.file("Dfam HMM", classification["dfam_db"], "classification.dfam_db")
            checks.executable("NHMMER", classification.get("nhmmer_bin", "nhmmer"), "classification.nhmmer_bin")
        return

    from drayte.classification.dfam_validation import options

    try:
        cfg = options(classification["dfam"])
    except (TypeError, ValueError) as error:
        checks.problem(str(error))
        return
    if not cfg["enabled"]:
        return

    warn = cfg["on_missing"] == "warn"
    checks.executable("NHMMER", cfg.get("nhmmer_bin", "nhmmer"), "classification.dfam.nhmmer_bin", warn=warn)
    database = checks.text(cfg.get("database"), "classification.dfam.database",
                           warn=warn and (cfg.get("database") is None or isinstance(cfg.get("database"), str)))
    if database is None:
        return
    path = Path(database).expanduser().resolve()
    if path.is_dir():
        checks.executable("FamDB", cfg.get("famdb_bin", "famdb.py"), "classification.dfam.famdb_bin", warn=warn)
        if not cfg.get("taxon"):
            checks.problem("classification.dfam.taxon: set the lineage to export from the FamDB directory.")
        try:
            partitions = sorted(path.glob("*.h5"))
            if not partitions:
                checks.problem(f"No uncompressed FamDB *.h5 partitions in {path} [classification.dfam.database].", warn=warn)
            for partition in partitions:
                checks.file("FamDB partition", str(partition), "classification.dfam.database", warn=warn)
        except OSError as error:
            checks.problem(f"Cannot read FamDB directory {path}: {error} [classification.dfam.database].", warn=warn)
    else:
        checks.file("Dfam HMM", str(path), "classification.dfam.database", warn=warn and not path.exists())
        if path.is_file():
            if path.suffix == ".h5":
                checks.problem("classification.dfam.database: use the FamDB directory, not an individual .h5 partition.")
            elif cfg.get("taxon"):
                checks.problem("classification.dfam.taxon: use null for an already exported HMM; taxonomic filtering needs a FamDB directory.")


def check_preflight(config: PipelineConfig) -> PreflightReport:
    """Check the configured workflow, without trusting or changing checkpoints.

    Commands are located, not launched. A resumed run is checked against the same
    tool requirements as a new run; stage completion is still validated later.
    """
    checks = _Checks()
    extra = config.extra
    checks.file("Genome FASTA", config.genome, "genome")
    try:
        if not config.outdir_path.is_absolute():
            checks.problem("outdir: use an absolute output directory path; extension changes its working directory.")
        destination = config.outdir_path.resolve()
        while not destination.exists():
            destination = destination.parent
        if not destination.is_dir() or not os.access(destination, os.W_OK | os.X_OK):
            checks.problem(f"Output directory or its existing parent is not writable: {destination} [outdir].")
    except OSError as error:
        checks.problem(f"Cannot access output directory {config.outdir}: {error} [outdir].")

    # Commands run by Python directly, independently of third-party configuration.
    for command in ("perl", "makeblastdb", "blastn", "faToTwoBit"):
        checks.executable(command, command, "PATH", changed_cwd=command == "perl")
    checks.executable("RepeatMasker", extra.get("repeatmasker_bin", "RepeatMasker"), "repeatmasker_bin", changed_cwd=True)
    checks.executable("RepeatClassifier", extra.get("repeatclassifier_bin", "RepeatClassifier"), "repeatclassifier_bin", changed_cwd=True)
    checks.executable("DIAMOND", extra.get("diamond_bin", "diamond"), "diamond_bin")

    modeler = checks.directory(extra.get("repeatmodeler_dir"), "repeatmodeler_dir", absolute=True)
    supplied_library = extra.get("repeatmodeler_library")
    # discovery.run indexes this setting even when an external library is used.
    scout = checks.text(extra.get("repeatscout_dir"), "repeatscout_dir")
    if supplied_library:
        library = checks.text(supplied_library, "repeatmodeler_library / --repeatmodeler-library")
        if library is not None:
            checks.file("RepeatModeler library", str(Path(library).expanduser()), "repeatmodeler_library / --repeatmodeler-library")
    else:
        if modeler is not None:
            for name in ("BuildDatabase", "RepeatModeler"):
                checks.executable(name, str(modeler / name), "repeatmodeler_dir")
        if scout is not None:
            for name in ("RepeatScout", "build_lmer_table"):
                checks.executable(name, str(Path(scout) / name), "repeatscout_dir")

    checks.file("Extension Perl script", extra.get("repeatmodeler_extend_script"), "repeatmodeler_extend_script", absolute=True)
    if modeler is not None:
        checks.file("Consensus extension helper", str(modeler / "util/alignAndCallConsensus.pl"), "repeatmodeler_dir")
        checks.executable("Linup", str(modeler / "util/Linup"), "repeatmodeler_dir")
        for matrix in ("crossmatch/18p41g.matrix", "ncbi/nt/25p41g-Hpad.matrix"):
            checks.file("Extension scoring matrix", str(modeler / "Matrices" / matrix), "repeatmodeler_dir")
        renderer = modeler / "util/visualizeAlignPNG.py"
        if renderer.is_file():
            checks.executable("Alignment renderer Python", "python", "PATH", warn=True, changed_cwd=True)
        elif not (modeler / "util/visualizeAlignPNG.pl").is_file():
            checks.problem(f"No alignment PNG renderer in {modeler / 'util'}; extension plots may be unavailable.", warn=True)

    # The supplied .2bit workflow does not use this directory directly. Linup
    # reads its UCSC location from RepeatModeler's own RepModelConfig instead.
    ucsc = checks.text(extra.get("ucsc_tools_dir"), "ucsc_tools_dir")
    if ucsc is not None:
        for name in ("faToTwoBit", "twoBitInfo", "twoBitToFa"):
            checks.executable(name, str(Path(ucsc) / name), "ucsc_tools_dir; also check RepeatModeler RepModelConfig", warn=True)

    for setting, names in (
        ("repeat_after_me_dir", ("RAMExtend",)),
        ("phrap_dir", ("cross_match",)),
    ):
        directory = checks.directory(extra.get(setting), setting, absolute=True)
        if directory is not None:
            for name in names:
                checks.executable(name, str(directory / name), setting)

    repeatpeps = checks.text(extra.get("repeatpeps_db_dir"), "repeatpeps_db_dir")
    if repeatpeps is not None:
        checks.file("RepeatPeps", str(Path(repeatpeps) / "RepeatPeps.lib"), "repeatpeps_db_dir")

    structure = checks.mapping(extra.get("structure", {}), "structure")
    heliano = checks.mapping(structure.get("heliano", {}), "structure.heliano")
    if heliano.get("enabled", False):
        checks.executable("HELIANO", heliano.get("bin", "heliano"), "structure.heliano.bin", changed_cwd=True)
        checks.executable("HELIANO MMseqs", heliano.get("mmseqs_bin", "mmseqs"), "structure.heliano.mmseqs_bin")
        for command in ("bedtools", "rnabob", "cd-hit-est", "mafft", "hmmsearch", "getorf", "gt", "dialign2-2", "Rscript", "bash"):
            checks.executable(command, command, "structure.heliano.enabled; PATH", changed_cwd=True)

    classification = checks.mapping(extra.get("classification", {}), "classification")
    # The integrated stage selects evidence by the configured references. Its
    # historical enabled/preliminary/final flags do not disable this stage.
    if classification.get("pfam_db"):
        checks.executable("HMMSCAN", classification.get("hmmscan_bin", "hmmscan"), "classification.hmmscan_bin")
        checks.file("Pfam HMM", classification["pfam_db"], "classification.pfam_db")
        for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
            checks.file("Pressed Pfam index (run hmmpress on the HMM to create it)", str(classification["pfam_db"]) + suffix, "classification.pfam_db")
    _check_dfam(checks, classification)
    if classification.get("mmseqs_rescue_db"):
        checks.executable("Classification MMseqs", classification.get("mmseqs_bin", "mmseqs"), "classification.mmseqs_bin")
        reference = checks.text(classification["mmseqs_rescue_db"], "classification.mmseqs_rescue_db / --mmseqs-rescue-db")
        if reference is not None:
            checks.file("MMseqs rescue FASTA", str(Path(reference).expanduser()), "classification.mmseqs_rescue_db / --mmseqs-rescue-db")

    return checks.finish()
