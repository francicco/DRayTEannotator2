from drayte.extension import consensus


class FakeProcess:
    stdout = ()

    def wait(self):
        return 0


def test_extension_command_passes_all_tool_directories(tmp_path, monkeypatch):
    captured = {}

    def fake_popen(cmd, **kwargs):
        captured["cmd"] = cmd
        captured["kwargs"] = kwargs
        return FakeProcess()

    monkeypatch.setattr(consensus.subprocess, "Popen", fake_popen)

    consensus.run_extend_consensus(
        extend_script="/tools/davidExtendConsRAM.pl",
        genome_2bit="/data/genome.2bit",
        family_fasta="/data/family.fa",
        workdir=tmp_path / "work",
        logfile=tmp_path / "extend.log",
        ucsc_tools_dir="/tools/ucsc",
        repeatmodeler_dir="/tools/repeatmodeler",
        repeat_after_me_dir="/tools/repeat-after-me",
        phrap_dir="/tools/phrap",
    )

    assert captured["cmd"] == [
        "perl",
        "/tools/davidExtendConsRAM.pl",
        "--ucscToolsDir",
        "/tools/ucsc",
        "--repeatModelerDir",
        "/tools/repeatmodeler",
        "--repeatAfterMeDir",
        "/tools/repeat-after-me",
        "--phrapDir",
        "/tools/phrap",
        "-genome",
        "/data/genome.2bit",
        "-family",
        "/data/family.fa",
        "-outdir",
        ".",
    ]
