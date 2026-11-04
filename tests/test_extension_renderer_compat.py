"""Exercise the actual extension renderer dispatch without a genome-scale run."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


@pytest.mark.skipif(not shutil.which('perl'), reason='Perl is required')
@pytest.mark.parametrize('custom_python', [False, True])
def test_extension_uses_available_repeatmodeler_renderer(tmp_path, custom_python):
    source = (Path(__file__).parents[1] / 'davidExtendConsRAM.pl').read_text()
    function = source.split('sub run_alignment_visualizer {', 1)[1].split('\nmy $matrixDir', 1)[0]
    harness = tmp_path / 'dispatch.pl'
    harness.write_text('use strict;\nmy $repeatModelerDir = shift;\nsub run_alignment_visualizer {' + function + '\nrun_alignment_visualizer("alignment with spaces.out", "genome with spaces.2bit");\n')
    tools = tmp_path / 'RepeatModeler with spaces' / 'util'
    tools.mkdir(parents=True)
    (tools / 'visualizeAlignPNG.pl').write_text('print "perl\\n"; print join("\\n", @ARGV);\n')
    env = os.environ.copy()
    if custom_python:
        (tools / 'visualizeAlignPNG.py').write_text('import sys\nprint("python")\nprint("\\n".join(sys.argv[1:]))\n')
        bindir = tmp_path / 'bin'
        bindir.mkdir()
        (bindir / 'python').symlink_to(sys.executable)
        env['PATH'] = str(bindir) + os.pathsep + env['PATH']
    result = subprocess.run(['perl', str(harness), str(tools.parent)], text=True, capture_output=True, check=True, env=env)
    assert result.stdout.splitlines() == ['python' if custom_python else 'perl', '-genome', 'genome with spaces.2bit', '-outfile', 'alignment with spaces.out']
    assert not result.stderr
