from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from barcode_kit.exceptions import PhylogenyError
from barcode_kit.phylogeny import (
    AlignmentProgram,
    SubprocessAlignmentRunner,
    SubprocessTreeRunner,
    SubprocessTreeShrinkRunner,
    SubprocessTrimalRunner,
    TreeShrinkResult,
    run_tree_shrink_qc,
)


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

    result = SubprocessAlignmentRunner().align(
        input_path,
        output_path,
        program=AlignmentProgram.MAFFT,
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


def test_alignment_runner_runs_muscle_with_output_argument(tmp_path: Path, monkeypatch):
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "aligned.fasta"
    input_path.write_text(">a\nACGT\n", encoding="utf-8")
    commands: list[list[str]] = []

    def fake_which(name: str) -> str:
        return f"/fake/{name}"

    def fake_run(command, **kwargs):
        commands.append(command)
        output_path.write_text(">a\nACGT\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("barcode_kit.phylogeny.shutil.which", fake_which)
    monkeypatch.setattr("barcode_kit.phylogeny.subprocess.run", fake_run)

    SubprocessAlignmentRunner().align(
        input_path,
        output_path,
        program=AlignmentProgram.MUSCLE,
        threads=2,
    )

    assert commands == [
        [
            "/fake/muscle",
            "-align",
            str(input_path),
            "-output",
            str(output_path),
            "-quiet",
            "-nt",
            "-threads",
            "2",
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

    result = SubprocessTrimalRunner().trim(input_path, output_path)

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

    result = SubprocessTreeRunner().build_tree(
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
        SubprocessAlignmentRunner().align(input_path, output_path)


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

    result = SubprocessTreeShrinkRunner().detect_outliers(
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
        SubprocessTreeShrinkRunner().detect_outliers(tree_path, tmp_path / "treeshrink")


def test_run_tree_shrink_qc_filters_original_fasta_after_alignment_and_tree(tmp_path: Path):
    input_path = tmp_path / "input.fasta"
    output_path = tmp_path / "filtered.fasta"
    workdir = tmp_path / "qc"
    input_path.write_text(">keep\nACGT\n>bad\nACGA\n", encoding="utf-8")

    class FakeAlignmentRunner:
        def align(self, input_path, output_path, *, program, threads):
            output_path.write_text(Path(input_path).read_text(encoding="utf-8"), encoding="utf-8")
            return output_path

    class FakeTreeRunner:
        def __init__(self):
            self.calls = []

        def build_tree(self, input_path, output_path, *, bootstrap):
            self.calls.append((input_path, output_path, bootstrap))
            output_path.write_text("(keep:0.1,bad:1.5);\n", encoding="utf-8")
            return output_path

    class FakeTreeShrinkRunner:
        def detect_outliers(
            self,
            tree_path,
            output_dir,
            *,
            output_prefix="output",
            quantile=0.05,
            max_removed=None,
        ):
            self.call = (tree_path, output_dir, output_prefix, quantile, max_removed)
            output_dir.mkdir(parents=True, exist_ok=True)
            removed_path = output_dir / f"{output_prefix}.txt"
            removed_path.write_text("bad\t\n", encoding="utf-8")
            return TreeShrinkResult(
                removed_taxa={"bad"},
                output_dir=output_dir,
                removed_taxa_path=removed_path,
            )

    tree_runner = FakeTreeRunner()
    tree_shrink_runner = FakeTreeShrinkRunner()
    result = run_tree_shrink_qc(
        input_path,
        output_path,
        workdir,
        bootstrap=1000,
        max_removed=6,
        alignment_runner=FakeAlignmentRunner(),
        tree_runner=tree_runner,
        tree_shrink_runner=tree_shrink_runner,
    )

    assert result.removed_taxa == {"bad"}
    assert result.alignment_path == workdir / "mafft.fasta"
    assert result.tree_path == workdir / "iqtree.tree"
    assert tree_runner.calls == [(workdir / "mafft.fasta", workdir / "iqtree.tree", 1000)]
    assert tree_shrink_runner.call == (
        workdir / "iqtree.tree",
        workdir / "treeshrink",
        "output",
        0.05,
        6,
    )
    assert output_path.read_text(encoding="utf-8") == ">keep\nACGT\n"
