import logging
import shutil
from concurrent.futures import ThreadPoolExecutor
from textwrap import dedent

import pytest

from drayte.extension.consensus import run_extend_consensus


@pytest.mark.skipif(shutil.which("perl") is None, reason="Perl is unavailable")
def test_parallel_extension_subprocesses_isolate_files_and_logs(tmp_path, caplog):
    """Exercise the real wrapper, with two children meeting at a bounded barrier."""
    barrier = tmp_path / "barrier"
    barrier.mkdir()
    script = tmp_path / "fake_extend.pl"
    script.write_text(
        dedent(
            r"""
            use strict;
            use warnings;
            use Cwd qw(getcwd);
            use File::Basename qw(basename dirname);
            use Getopt::Long qw(GetOptions);
            use Time::HiRes qw(time sleep);

            my %options;
            GetOptions(\%options, 'ucscToolsDir=s', 'repeatModelerDir=s',
                'repeatAfterMeDir=s', 'phrapDir=s', 'genome=s',
                'family=s', 'outdir=s') or die "bad arguments";
            die "unexpected output directory" unless $options{outdir} eq '.';
            my $name = basename($options{family});
            $name =~ s/\.fa$//;
            $| = 1;

            # Both jobs use the same intermediate names in different cwd values.
            open my $scratch, '>', 'scratch.txt' or die $!;
            print {$scratch} "$name\n";
            close $scratch;
            my $barrier = dirname($0) . '/barrier';
            open my $ready, '>', "$barrier/$name.ready" or die $!;
            close $ready;
            my $deadline = time() + 5;
            while (1) {
                my @ready = glob("$barrier/*.ready");
                last if @ready == 2;
                die "children did not overlap" if time() > $deadline;
                sleep(0.01);
            }

            open my $check, '<', 'scratch.txt' or die $!;
            my $owner = <$check>;
            close $check;
            die "shared scratch file was overwritten" unless $owner eq "$name\n";
            open my $rep, '>', 'rep' or die $!;
            print {$rep} ">$name\nACGT\n";
            close $rep;
            print "family=$name cwd=" . getcwd() . "\n";
            print STDERR "stderr family=$name\n";
            """
        )
    )
    genome = tmp_path / "genome.2bit"
    genome.write_text("shared read-only genome\n")
    families = ["family-A#Unknown", "family-B#DNA_TIR"]
    paths = {}
    for family in families:
        family_fasta = tmp_path / f"{family}.fa"
        family_fasta.write_text(f">{family}\nACGT\n")
        paths[family] = {
            "family_fasta": family_fasta,
            "workdir": tmp_path / "work" / family,
            "logfile": tmp_path / f"{family}.log",
        }

    caplog.set_level(logging.INFO, logger="drayte.extension.consensus")
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(
                run_extend_consensus,
                extend_script=script,
                genome_2bit=genome,
                ucsc_tools_dir=tmp_path,
                repeatmodeler_dir=tmp_path,
                repeat_after_me_dir=tmp_path,
                phrap_dir=tmp_path,
                **paths[family],
            )
            for family in families
        ]
        for future in futures:
            future.result(timeout=10)

    messages = [record.getMessage() for record in caplog.records]
    for family in families:
        workdir = paths[family]["workdir"]
        assert (workdir / "scratch.txt").read_text() == f"{family}\n"
        assert (workdir / "rep").read_text() == f">{family}\nACGT\n"
        log = paths[family]["logfile"].read_text()
        assert f"family={family} cwd={workdir}" in log
        assert f"stderr family={family}" in log
        other_family = next(name for name in families if name != family)
        assert other_family not in log
        assert f"[extend:{family}] family={family} cwd={workdir}" in messages
        assert f"[extend:{family}] stderr family={family}" in messages

    assert not (tmp_path / "scratch.txt").exists()
    assert not (tmp_path / "rep").exists()
    assert genome.read_text() == "shared read-only genome\n"
