from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Protocol

from Bio import SeqIO

from barcode_kit.exceptions import PhylogenyError


class AlignmentProgram(StrEnum):
    MAFFT = "mafft"
    MUSCLE = "muscle"
    CLUSTALO = "clustalo"


class TreeProgram(StrEnum):
    FASTTREE = "fasttree"
    VERYFASTTREE = "veryfasttree"
    IQTREE = "iqtree"
    RAXMLNG = "raxmlng"


ALIGNMENT_COMMANDS = {
    AlignmentProgram.MAFFT: "mafft",
    AlignmentProgram.MUSCLE: "muscle",
    AlignmentProgram.CLUSTALO: "clustalo",
}
TREE_COMMANDS = {
    TreeProgram.FASTTREE: "FastTree",
    TreeProgram.VERYFASTTREE: "VeryFastTree",
    TreeProgram.IQTREE: "iqtree3",
    TreeProgram.RAXMLNG: "raxml-ng",
}
TRIMAL_COMMAND = "trimal"
TREESHRINK_COMMAND = "run_treeshrink.py"


@dataclass(frozen=True)
class TreeShrinkResult:
    removed_taxa: set[str]
    output_dir: Path
    removed_taxa_path: Path


@dataclass(frozen=True)
class TreeShrinkQcConfig:
    threads: int = 1
    quantile: float = 0.05


@dataclass(frozen=True)
class TreeShrinkQcResult:
    removed_taxa: set[str]
    output_fasta: Path
    alignment_path: Path
    tree_path: Path
    tree_shrink_output_dir: Path
    removed_taxa_path: Path


class AlignmentRunner(Protocol):
    def align(
        self,
        input_path: Path,
        output_path: Path,
        *,
        program: AlignmentProgram,
        threads: int,
    ) -> Path:
        ...


class TreeRunner(Protocol):
    def build_tree(
        self,
        input_path: Path,
        output_path: Path,
        *,
        program: TreeProgram,
        threads: int,
    ) -> Path:
        ...


class TreeShrinkRunner(Protocol):
    def detect_outliers(
        self,
        tree_path: Path,
        output_dir: Path,
        *,
        output_prefix: str = "output",
        quantile: float = 0.05,
    ) -> TreeShrinkResult:
        ...


class SubprocessAlignmentRunner:
    def align(
        self,
        input_path: Path,
        output_path: Path,
        *,
        program: AlignmentProgram = AlignmentProgram.MAFFT,
        threads: int = 1,
    ) -> Path:
        executable = _find_executable(ALIGNMENT_COMMANDS[program])
        input_path = Path(input_path)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        threads_text = str(max(1, threads))

        if program is AlignmentProgram.MAFFT:
            command = [
                executable,
                "--auto",
                "--quiet",
                "--nuc",
                "--thread",
                threads_text,
                str(input_path),
            ]
            with output_path.open("w", encoding="utf-8") as output_handle:
                _run_command(command, stdout=output_handle)
        elif program is AlignmentProgram.MUSCLE:
            command = [
                executable,
                "-align",
                str(input_path),
                "-output",
                str(output_path),
                "-quiet",
                "-nt",
                "-threads",
                threads_text,
            ]
            _run_command(command)
        else:
            command = [
                executable,
                "-i",
                str(input_path),
                "-o",
                str(output_path),
                "--auto",
                "--force",
                "--seqtype=DNA",
                f"--threads={threads_text}",
            ]
            _run_command(command)

        return _ensure_output(output_path, "multiple sequence alignment")


class SubprocessTrimalRunner:
    def trim(
        self,
        input_path: Path,
        output_path: Path,
        *,
        automated1: bool = True,
    ) -> Path:
        executable = _find_executable(TRIMAL_COMMAND)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        command = [
            executable,
            "-in",
            str(input_path),
            "-out",
            str(output_path),
        ]
        if automated1:
            command.append("-automated1")
        _run_command(command)
        return _ensure_output(output_path, "trimAl sequence trimming")


class SubprocessTreeRunner:
    def build_tree(
        self,
        input_path: Path,
        output_path: Path,
        *,
        program: TreeProgram = TreeProgram.FASTTREE,
        bootstrap: int = 0,
        threads: int = 1,
    ) -> Path:
        executable = _find_executable(TREE_COMMANDS[program])
        input_path = Path(input_path)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        tool_output_path = _tool_tree_output_path(input_path, program)
        command = _tree_command(
            executable,
            input_path,
            tool_output_path,
            program,
            bootstrap,
            max(1, threads),
        )
        _run_command(command)
        _ensure_output(tool_output_path, "phylogenetic tree reconstruction")
        shutil.copyfile(tool_output_path, output_path)
        return _ensure_output(output_path, "phylogenetic tree reconstruction")


class SubprocessTreeShrinkRunner:
    def detect_outliers(
        self,
        tree_path: Path,
        output_dir: Path,
        *,
        output_prefix: str = "output",
        quantile: float = 0.05,
    ) -> TreeShrinkResult:
        executable = _find_executable(TREESHRINK_COMMAND)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        command = [
            executable,
            "-t",
            str(tree_path),
            "-o",
            str(output_dir),
            "-O",
            output_prefix,
            "-q",
            str(quantile),
            "-f",
        ]
        _run_command(command)
        removed_taxa_path = _ensure_output(output_dir / f"{output_prefix}.txt", "TreeShrink")
        return TreeShrinkResult(
            removed_taxa=_read_removed_taxa(removed_taxa_path),
            output_dir=output_dir,
            removed_taxa_path=removed_taxa_path,
        )


def run_tree_shrink_qc(
    input_fasta: Path,
    output_fasta: Path,
    workdir: Path,
    *,
    threads: int = 1,
    quantile: float = 0.05,
    alignment_runner: AlignmentRunner | None = None,
    tree_runner: TreeRunner | None = None,
    tree_shrink_runner: TreeShrinkRunner | None = None,
) -> TreeShrinkQcResult:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    alignment_path = workdir / "mafft.fasta"
    tree_path = workdir / "iqtree.tree"
    tree_shrink_output_dir = workdir / "treeshrink"
    alignment_runner = alignment_runner or SubprocessAlignmentRunner()
    tree_runner = tree_runner or SubprocessTreeRunner()
    tree_shrink_runner = tree_shrink_runner or SubprocessTreeShrinkRunner()

    alignment_runner.align(
        input_fasta,
        alignment_path,
        program=AlignmentProgram.MAFFT,
        threads=threads,
    )
    tree_runner.build_tree(
        alignment_path,
        tree_path,
        program=TreeProgram.IQTREE,
        threads=threads,
    )
    tree_shrink_result = tree_shrink_runner.detect_outliers(
        tree_path,
        tree_shrink_output_dir,
        quantile=quantile,
    )
    _write_filtered_fasta(input_fasta, output_fasta, tree_shrink_result.removed_taxa)
    return TreeShrinkQcResult(
        removed_taxa=tree_shrink_result.removed_taxa,
        output_fasta=output_fasta,
        alignment_path=alignment_path,
        tree_path=tree_path,
        tree_shrink_output_dir=tree_shrink_output_dir,
        removed_taxa_path=tree_shrink_result.removed_taxa_path,
    )


def _tree_command(
    executable: str,
    input_path: Path,
    tool_output_path: Path,
    program: TreeProgram,
    bootstrap: int,
    threads: int,
) -> list[str]:
    if program is TreeProgram.IQTREE:
        command = [executable, "-s", str(input_path), "-redo"]
        if bootstrap:
            command.extend(["-B", str(bootstrap)])
        if threads > 1:
            command.extend(["-T", "AUTO", "-ntmax", str(threads)])
        else:
            command.extend(["-T", "1"])
        return command

    if program is TreeProgram.RAXMLNG:
        command = [
            executable,
            "--msa",
            str(input_path),
            "--msa-format",
            "FASTA",
            "--model",
            "GTR+G",
            "--redo",
        ]
        if bootstrap:
            command.extend(["--all", "--bs-trees", str(bootstrap)])
        else:
            command.append("--search")
        if threads > 1:
            command.extend(["--threads", f"auto{{{threads}}}", "--workers", "auto"])
        else:
            command.extend(["--threads", "1"])
        return command

    command = [executable, "-out", str(tool_output_path), "-gtr"]
    if bootstrap:
        command.extend(["-boot", str(bootstrap)])
    else:
        command.append("-nosupport")
    if program is TreeProgram.VERYFASTTREE and threads > 1:
        command.extend(["-threads", str(threads)])
    command.extend(["-nt", str(input_path)])
    return command


def _tool_tree_output_path(input_path: Path, program: TreeProgram) -> Path:
    if program is TreeProgram.IQTREE:
        return Path(str(input_path) + ".treefile")
    if program is TreeProgram.RAXMLNG:
        return Path(str(input_path) + ".raxml.bestTree")
    if program is TreeProgram.VERYFASTTREE:
        return Path(str(input_path) + ".veryfasttree.tre")
    return Path(str(input_path) + ".fasttree.tre")


def _find_executable(command: str) -> str:
    executable = shutil.which(command)
    if executable is None:
        raise PhylogenyError(f"Unable to find {command} executable")
    return executable


def _run_command(command: list[str], **kwargs) -> None:
    try:
        completed = subprocess.run(
            command,
            check=False,
            stdout=kwargs.pop("stdout", subprocess.DEVNULL),
            stderr=subprocess.DEVNULL,
            **kwargs,
        )
    except FileNotFoundError as error:
        raise PhylogenyError(f"Unable to run {command[0]}") from error
    if completed.returncode != 0:
        raise PhylogenyError(f"{Path(command[0]).name} failed with exit code {completed.returncode}")


def _ensure_output(path: Path, operation: str) -> Path:
    if not path.exists():
        raise PhylogenyError(f"{operation} did not create {path}")
    return path


def _read_removed_taxa(path: Path) -> set[str]:
    taxa: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        taxa.update(token for token in line.split() if token)
    return taxa


def _write_filtered_fasta(input_fasta: Path, output_fasta: Path, removed_taxa: set[str]) -> None:
    records = [record for record in SeqIO.parse(str(input_fasta), "fasta") if record.id not in removed_taxa]
    with output_fasta.open("w", encoding="utf-8") as handle:
        SeqIO.write(records, handle, "fasta")
