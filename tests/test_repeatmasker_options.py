import json
import logging
import shutil
import subprocess
from types import SimpleNamespace

import pytest

from drayte import step1_repmodannotation as step1
from drayte.pipeline import final_annotation
from drayte.utils.repeatmasker_checkpoint import checkpoint_path


# Use RepeatMasker's actual Perl option semantics: gc=s consumes the next
# argument, including "-gff". Checking only that "-gff" occurs in argv misses it.
PERL_PARSE_OPTIONS = r"""
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use JSON::PP qw(encode_json);
Getopt::Long::Configure("noignorecase", "bundling_override");
my @argv = @ARGV;
my %options;
GetOptions(\%options, qw(gc=s gccalc gff xsmall norna a pa=i lib=s s))
    or die "Invalid RepeatMasker options\n";
my $parsed = {argv => \@argv, options => \%options, inputs => \@ARGV};
"""


@pytest.fixture(params=["discovery", "final"])
def stage(tmp_path, request):
    perl = shutil.which("perl")
    if perl is None:
        pytest.skip("Perl is needed to exercise RepeatMasker's option parser")
    script = tmp_path / "RepeatMasker"
    script.write_text(
        f"#!{perl}\n"
        + PERL_PARSE_OPTIONS
        + r"""
print "PARSED=" . encode_json($parsed) . "\n";
@ARGV == 1 or die "Expected one genome\n";
my @suffixes = qw(.out .tbl .masked .cat.gz .align);
push @suffixes, ".gff" if $options{gff};
for my $suffix (@suffixes) {
    open my $handle, ">", $ARGV[0] . $suffix or die $!;
    print {$handle} "RepeatMasker output $suffix\n";
    close $handle or die $!;
}
"""
    )
    script.chmod(0o755)
    root = tmp_path / "run"
    genome = root / "discovery" / "assemblies_dir" / "T.fa"
    genome.parent.mkdir(parents=True)
    genome.write_text(">chr1\nACGTACGT\n")
    library = tmp_path / "library.fa"
    library.write_text(">family#RC/Helitron\nACGT\n")
    logger = logging.getLogger("test.repeatmasker.options")
    if request.param == "discovery":
        outdir = root / "discovery" / "rmasker_dir"

        def run():
            return step1.run_step1(
                genome=genome,
                outdir=root / "discovery",
                species="T",
                threads=3,
                repeatmodeler_dir=tmp_path / "unused-repeatmodeler",
                repeatscout_dir=tmp_path / "unused-repeatscout",
                repeatmasker_bin=str(script),
                repeatmodeler_library=library,
                logger=logger,
            )
    else:
        outdir = root / "final_annotation"
        config = SimpleNamespace(
            outdir_path=root,
            species="T",
            threads=3,
            extra={"repeatmasker_bin": str(script)},
        )

        def run():
            return final_annotation.run(
                config, {"classified_library": str(library)}, logger
            )

    return SimpleNamespace(
        run=run, outdir=outdir, perl=perl, name=request.param, genome=genome
    )


def parsed_attempts(stage):
    log = (stage.outdir / "T.RepeatMasker.log").read_text()
    return [
        json.loads(line.removeprefix("PARSED="))
        for line in log.splitlines()
        if line.startswith("PARSED=")
    ]


def test_stage_gff_and_default_gc_policy_match_checkpoint(stage):
    stage.run()
    parsed = parsed_attempts(stage)[0]
    actual_options = parsed["options"]
    assert actual_options["gff"] == 1
    assert "gc" not in actual_options
    assert "gccalc" not in actual_options
    assert all(actual_options[flag] == 1 for flag in ("s", "xsmall", "norna"))
    assert actual_options.get("a", 0) == (stage.name == "final")
    assert actual_options["pa"] == 3
    assert parsed["inputs"] == [str(stage.genome)]
    assert (stage.outdir / "T.fa.gff").is_file()

    checkpoint = json.loads(checkpoint_path(stage.outdir).read_text())
    checkpoint_argv = checkpoint["input_signature"]["options"][1:]
    completed = subprocess.run(
        [stage.perl, "-e", PERL_PARSE_OPTIONS + "print encode_json($parsed);",
         "--", *checkpoint_argv],
        check=True,
        capture_output=True,
        text=True,
    )
    checkpoint_options = json.loads(completed.stdout)["options"]
    assert checkpoint_options == {
        key: value for key, value in actual_options.items() if key not in ("lib", "pa")
    }

    # The matching signature must still permit normal resume.
    stage.run()
    assert len(parsed_attempts(stage)) == 1


def test_checkpoint_with_old_bare_gc_option_is_rerun(stage):
    stage.run()
    path = checkpoint_path(stage.outdir)
    old_checkpoint = json.loads(path.read_text())
    options = old_checkpoint["input_signature"]["options"]
    options.insert(options.index("-gff"), "-gc")
    path.write_text(json.dumps(old_checkpoint))

    stage.run()

    assert len(parsed_attempts(stage)) == 2
    assert "-gc" not in json.loads(path.read_text())["input_signature"]["options"]
