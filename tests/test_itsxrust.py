from __future__ import annotations

import subprocess
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from barcode_kit.config import ItsxrustConfig
from barcode_kit.itsxrust import ItsxrustInput, SubprocessItsxrustRunner
from barcode_kit.models import Marker


def test_subprocess_itsxrust_runner_extract_many_uses_one_process_and_maps_missing_outputs(
    tmp_path: Path,
    monkeypatch,
):
    hmm_path = tmp_path / "T.hmm"
    hmm_path.write_text("hmm", encoding="utf-8")
    commands: list[list[str]] = []

    def fake_which(name: str) -> str:
        return f"/fake/{name}"

    def fake_run(command, check, capture_output, text):
        commands.append(command)
        output_path = Path(command[command.index("--output") + 1])
        output_path.write_text(">ITS000001.1\nAAAACCCCGGGG\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("barcode_kit.itsxrust.shutil.which", fake_which)
    monkeypatch.setattr("barcode_kit.itsxrust.subprocess.run", fake_run)

    runner = SubprocessItsxrustRunner()
    results = runner.extract_many(
        [
            ItsxrustInput("ITS000001.1", SeqRecord(Seq("AAAACCCCGGGG"))),
            ItsxrustInput("ITS000002.1", SeqRecord(Seq("TTTTCCCCAAAA"))),
        ],
        Marker.ITS,
        hmm_path,
    )

    assert len(commands) == 1
    assert results["ITS000001.1"].sequence == Seq("AAAACCCCGGGG")
    assert results["ITS000001.1"].fallback_reason is None
    assert results["ITS000002.1"].sequence is None
    assert results["ITS000002.1"].fallback_reason == "no_anchor_full"
    command = commands[0]
    assert command[command.index("--inc-e") + 1] == "0.01"
    assert command[command.index("--min-anchor-score") + 1] == "8"
    assert command[command.index("--max-per-anchor") + 1] == "20"
    assert command[command.index("--max-anchor-evalue") + 1] == "0.01"


def test_subprocess_itsxrust_runner_uses_configured_anchor_params(
    tmp_path: Path,
    monkeypatch,
):
    hmm_path = tmp_path / "T.hmm"
    hmm_path.write_text("hmm", encoding="utf-8")
    commands: list[list[str]] = []

    def fake_which(name: str) -> str:
        return f"/fake/{name}"

    def fake_run(command, check, capture_output, text):
        commands.append(command)
        output_path = Path(command[command.index("--output") + 1])
        output_path.write_text(">ITS000001.1\nAAAACCCCGGGG\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("barcode_kit.itsxrust.shutil.which", fake_which)
    monkeypatch.setattr("barcode_kit.itsxrust.subprocess.run", fake_run)

    runner = SubprocessItsxrustRunner(
        config=ItsxrustConfig(
            inc_e=0.02,
            min_anchor_score=10,
            max_per_anchor=25,
            max_anchor_evalue=0.02,
        )
    )
    runner.extract_many(
        [ItsxrustInput("ITS000001.1", SeqRecord(Seq("AAAACCCCGGGG")))],
        Marker.ITS,
        hmm_path,
    )

    command = commands[0]
    assert command[command.index("--inc-e") + 1] == "0.02"
    assert command[command.index("--min-anchor-score") + 1] == "10"
    assert command[command.index("--max-per-anchor") + 1] == "25"
    assert command[command.index("--max-anchor-evalue") + 1] == "0.02"
