"""Execute the extension script's real legacy/upstream visualizer dispatch."""

import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


PERL = shutil.which("perl")
SCRIPT = Path(__file__).resolve().parents[1] / "davidExtendConsRAM.pl"
pytestmark = pytest.mark.skipif(PERL is None, reason="Perl is unavailable")


def render(tmp_path, modeler, *, alignment="out", extra_env=None):
    environment = os.environ.copy()
    environment.update(extra_env or {})
    genome = tmp_path / "genome with spaces;literal.2bit"
    return subprocess.run(
        [PERL, "-e", "my $script = shift @ARGV; require $script; render_alignment_png(@ARGV);",
         str(SCRIPT), str(modeler), str(genome), alignment],
        cwd=tmp_path, env=environment, capture_output=True, text=True,
    )


def add_perl_renderer(modeler, *, exit_code=0):
    utilities = modeler / "util"
    utilities.mkdir(parents=True, exist_ok=True)
    renderer = utilities / "visualizeAlignPNG.pl"
    renderer.write_text(
        "use strict;\n"
        "open my $args, '>', 'perl.args' or die $!;\n"
        "print {$args} join(qq{\\n}, @ARGV); close $args;\n"
        "open my $png, '>', 'img.png' or die $!; print {$png} 'test PNG'; close $png;\n"
        f"exit {exit_code};\n"
    )
    return renderer


def test_loading_extension_script_does_not_start_extension(tmp_path):
    result = subprocess.run(
        [PERL, "-e", "my $script = shift @ARGV; require $script;", str(SCRIPT)],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == result.stderr == ""
    assert list(tmp_path.iterdir()) == []


def test_direct_extension_entrypoint_still_displays_help(tmp_path):
    result = subprocess.run(
        [PERL, str(SCRIPT), "--ucscToolsDir", str(tmp_path),
         "--repeatModelerDir", str(tmp_path), "--repeatAfterMeDir", str(tmp_path),
         "--phrapDir", str(tmp_path), "-h"],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "davidExtendConsRAM.pl -genome" in result.stdout
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("alignment", ["out", "repam-newrep.out"])
def test_upstream_perl_visualizer_receives_exact_arguments_and_writes_image(tmp_path, alignment):
    modeler = tmp_path / "RepeatModeler with spaces"
    add_perl_renderer(modeler)
    result = render(tmp_path, modeler, alignment=alignment)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "perl.args").read_text().splitlines() == [
        "-genome", str(tmp_path / "genome with spaces;literal.2bit"),
        "-outfile", alignment,
    ]
    assert (tmp_path / "img.png").read_bytes() == b"test PNG"


def test_legacy_python_visualizer_is_preferred_when_both_are_present(tmp_path):
    modeler = tmp_path / "RepeatModeler"
    add_perl_renderer(modeler)
    (modeler / "util/visualizeAlignPNG.py").write_text(
        "from pathlib import Path\nimport sys\n"
        "Path('python.args').write_text('\\n'.join(sys.argv[1:]))\n"
        "Path('img.png').write_bytes(b'legacy PNG')\n"
    )
    executables = tmp_path / "bin"
    executables.mkdir()
    (executables / "python").symlink_to(sys.executable)
    result = render(tmp_path, modeler, extra_env={
        "PATH": str(executables) + os.pathsep + os.environ.get("PATH", ""),
    })
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "python.args").read_text().splitlines() == [
        "-genome", str(tmp_path / "genome with spaces;literal.2bit"),
        "-outfile", "out",
    ]
    assert (tmp_path / "img.png").read_bytes() == b"legacy PNG"
    assert not (tmp_path / "perl.args").exists()


def test_missing_visualizer_has_an_actionable_error(tmp_path):
    result = render(tmp_path, tmp_path / "RepeatModeler")
    assert result.returncode != 0
    assert "No RepeatModeler alignment visualizer found" in result.stderr
    assert "visualizeAlignPNG.py or visualizeAlignPNG.pl" in result.stderr


def test_failed_visualizer_does_not_report_success(tmp_path):
    modeler = tmp_path / "RepeatModeler"
    add_perl_renderer(modeler, exit_code=7)
    result = render(tmp_path, modeler)
    assert result.returncode != 0
    assert "failed with exit code 7" in result.stderr

