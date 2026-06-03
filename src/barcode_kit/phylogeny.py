from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from barcode_kit.exceptions import PhylogenyError


__all__ = [
    "AlignmentRunner",
    "TreeRunner",
    "TreeShrinkRunner",
    "TreeShrinkResult",
    "TrimalRunner",
]


MAFFT_COMMAND = "mafft"
IQTREE_COMMAND = "iqtree3"
TRIMAL_COMMAND = "trimal"
TREESHRINK_COMMAND = "run_treeshrink.py"


@dataclass(frozen=True)
class TreeShrinkResult:
    removed_taxa: set[str]
    output_dir: Path
    removed_taxa_path: Path


class AlignmentRunner:
    def align(
        self,
        input_path: Path,
        output_path: Path,
        *,
        threads: int = 1,
    ) -> Path:
        executable = _find_executable(MAFFT_COMMAND)
        input_path = Path(input_path)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        threads_text = str(max(1, threads))
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

        return _ensure_output(output_path, "multiple sequence alignment")


class TrimalRunner:
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


class TreeRunner:
    def build_tree(
        self,
        input_path: Path,
        output_path: Path,
        *,
        bootstrap: int = 0,
    ) -> Path:
        executable = _find_executable(IQTREE_COMMAND)
        input_path = Path(input_path)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        tool_output_path = Path(str(input_path) + ".treefile")
        command = [executable, "-s", str(input_path), "-redo", "-m", "MFP"]
        if bootstrap:
            command.extend(["-B", str(bootstrap)])
        command.extend(["-T", "AUTO"])
        _run_command(command)
        _ensure_output(tool_output_path, "phylogenetic tree reconstruction")
        shutil.copyfile(tool_output_path, output_path)
        return _ensure_output(output_path, "phylogenetic tree reconstruction")


class TreeShrinkRunner:
    def detect_outliers(
        self,
        tree_path: Path,
        output_dir: Path,
        *,
        output_prefix: str = "output",
        quantile: float = 0.05,
        max_removed: int | None = None,
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
        ]
        if max_removed is not None:
            command.extend(["-k", str(max_removed)])
        command.extend(["-m", "per-gene", "-f"])
        _run_command(command)
        removed_taxa_path = _ensure_output(output_dir / f"{output_prefix}.txt", "TreeShrink")
        removed_taxa: set[str] = set()
        for line in removed_taxa_path.read_text(encoding="utf-8").splitlines():
            removed_taxa.update(token for token in line.split() if token)
        return TreeShrinkResult(
            removed_taxa=removed_taxa,
            output_dir=output_dir,
            removed_taxa_path=removed_taxa_path,
        )


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
