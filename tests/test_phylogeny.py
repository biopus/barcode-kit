from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

import barcode_kit.phylogeny as phylogeny_module
from barcode_kit.exceptions import PhylogenyError
from barcode_kit.phylogeny import (
    AlignmentRunner,
    TreeRunner,
    TreeShrinkRunner,
    TrimalRunner,
)


def test_phylogeny_module_only_exposes_tool_runners_and_tool_results():
    assert "AlignmentRunner" in phylogeny_module.__all__
    assert "AlignmentProgram" not in phylogeny_module.__all__
    assert "TreeRunner" in phylogeny_module.__all__
    assert "TreeShrinkRunner" in phylogeny_module.__all__
    assert "TrimalRunner" in phylogeny_module.__all__
    assert "run_tree_shrink_qc" not in phylogeny_module.__all__
    assert "TreeShrinkQcResult" not in phylogeny_module.__all__
    assert "SubprocessAlignmentRunner" not in phylogeny_module.__all__
    assert "SubprocessTreeRunner" not in phylogeny_module.__all__
    assert "SubprocessTreeShrinkRunner" not in phylogeny_module.__all__
    assert "SubprocessTrimalRunner" not in phylogeny_module.__all__
    assert not hasattr(phylogeny_module, "run_tree_shrink_qc")
    assert not hasattr(phylogeny_module, "TreeShrinkQcResult")
    assert not hasattr(phylogeny_module, "SubprocessAlignmentRunner")
    assert not hasattr(phylogeny_module, "SubprocessTreeRunner")
    assert not hasattr(phylogeny_module, "SubprocessTreeShrinkRunner")
    assert not hasattr(phylogeny_module, "SubprocessTrimalRunner")


def test_alignment_commands_only_include_mafft():
    assert phylogeny_module.MAFFT_COMMAND == "mafft"
    assert not hasattr(phylogeny_module, "AlignmentProgram")
    assert not hasattr(phylogeny_module, "ALIGNMENT_COMMANDS")


def test_alignment_runner_runs_mafft_and_writes_stdout_to_output(tmp_path: Path, monkeypatch):
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "aligned.fasta"
    input_path.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
    commands: list[list[str]] = []

    def fake_which(name: str) -> str:
        return f"/fake/{name}"

    def fake_run(command, **kwargs):
        commands.append(command)
        kwargs["stdout"].write(">a\nACGT\n>b\nACGA\n")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("barcode_kit.phylogeny.shutil.which", fake_which)
    monkeypatch.setattr("barcode_kit.phylogeny.subprocess.run", fake_run)

    result = AlignmentRunner().align(
        input_path,
        output_path,
        threads=4,
    )

    assert result == output_path
    assert output_path.read_text(encoding="utf-8") == ">a\nACGT\n>b\nACGA\n"
    assert commands == [
        [
            "/fake/mafft",
            "--auto",
            "--quiet",
            "--nuc",
            "--thread",
            "4",
            str(input_path),
        ]
    ]


def test_trimal_runner_runs_automated1(tmp_path: Path, monkeypatch):
    input_path = tmp_path / "aligned.fasta"
    output_path = tmp_path / "trimmed.fasta"
    input_path.write_text(">a\nACGT\n", encoding="utf-8")
    commands: list[list[str]] = []

    def fake_which(name: str) -> str:
        return f"/fake/{name}"

    def fake_run(command, **kwargs):
        commands.append(command)
        output_path.write_text(">a\nAC\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("barcode_kit.phylogeny.shutil.which", fake_which)
    monkeypatch.setattr("barcode_kit.phylogeny.subprocess.run", fake_run)

    result = TrimalRunner().trim(input_path, output_path)

    assert result == output_path
    assert commands == [
        ["/fake/trimal", "-in", str(input_path), "-out", str(output_path), "-automated1"]
    ]


def test_tree_runner_runs_iqtree_with_mfp_auto_threads_and_bootstrap(
    tmp_path: Path,
    monkeypatch,
):
    input_path = tmp_path / "trimmed.fasta"
    output_path = tmp_path / "tree.nwk"
    input_path.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
    commands: list[list[str]] = []

    def fake_which(name: str) -> str:
        return f"/fake/{name}"

    def fake_run(command, **kwargs):
        commands.append(command)
        (tmp_path / "trimmed.fasta.treefile").write_text("(a,b);\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("barcode_kit.phylogeny.shutil.which", fake_which)
    monkeypatch.setattr("barcode_kit.phylogeny.subprocess.run", fake_run)

    result = TreeRunner().build_tree(
        input_path,
        output_path,
        bootstrap=1000,
    )

    assert result == output_path
    assert output_path.read_text(encoding="utf-8") == "(a,b);\n"
    assert commands == [
        [
            "/fake/iqtree3",
            "-s",
            str(input_path),
            "-redo",
            "-m",
            "MFP",
            "-B",
            "1000",
            "-T",
            "AUTO",
        ]
    ]


def test_phylogeny_runners_raise_when_tool_is_unavailable(tmp_path: Path, monkeypatch):
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "aligned.fasta"
    input_path.write_text(">a\nACGT\n", encoding="utf-8")

    monkeypatch.setattr("barcode_kit.phylogeny.shutil.which", lambda name: None)

    with pytest.raises(PhylogenyError, match="Unable to find mafft executable"):
        AlignmentRunner().align(input_path, output_path)


def test_treeshrink_runner_runs_installed_command_and_reads_removed_taxa(
    tmp_path: Path,
    monkeypatch,
):
    tree_path = tmp_path / "input.tree"
    output_dir = tmp_path / "treeshrink"
    tree_path.write_text("(keep:0.1,bad:1.5);\n", encoding="utf-8")
    commands: list[list[str]] = []

    def fake_which(name: str) -> str:
        return f"/fake/{name}"

    def fake_run(command, **kwargs):
        commands.append(command)
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "output.txt").write_text("bad\t\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("barcode_kit.phylogeny.shutil.which", fake_which)
    monkeypatch.setattr("barcode_kit.phylogeny.subprocess.run", fake_run)

    result = TreeShrinkRunner().detect_outliers(
        tree_path,
        output_dir,
        quantile=0.01,
        max_removed=6,
    )

    assert result.removed_taxa == {"bad"}
    assert result.removed_taxa_path == output_dir / "output.txt"
    assert commands == [
        [
            "/fake/run_treeshrink.py",
            "-t",
            str(tree_path),
            "-o",
            str(output_dir),
            "-O",
            "output",
            "-q",
            "0.01",
            "-k",
            "6",
            "-m",
            "per-gene",
            "-f",
        ]
    ]


def test_treeshrink_runner_requires_installed_command(tmp_path: Path, monkeypatch):
    tree_path = tmp_path / "input.tree"
    tree_path.write_text("(a:0.1,b:0.2);\n", encoding="utf-8")
    monkeypatch.setattr("barcode_kit.phylogeny.shutil.which", lambda name: None)

    with pytest.raises(PhylogenyError, match="Unable to find run_treeshrink.py executable"):
        TreeShrinkRunner().detect_outliers(tree_path, tmp_path / "treeshrink")
