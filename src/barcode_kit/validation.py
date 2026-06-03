from __future__ import annotations

import json
from dataclasses import asdict, replace
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from barcode_kit.config import TreeShrinkConfig
from barcode_kit.models import BuildReportEntry, Marker, SequenceQuality
from barcode_kit.phylogeny import (
    AlignmentRunner,
    TreeRunner,
    TreeShrinkRunner,
)


__all__ = [
    "sequence_quality",
    "tree_shrink_qc",
]


def sequence_quality(sequence: Seq | str) -> SequenceQuality:
    seq = str(sequence).upper().replace(" ", "").replace("\n", "")
    length = len(seq)
    canonical_count = sum(seq.count(base) for base in "ACGT")
    gc_content = ((seq.count("G") + seq.count("C")) / canonical_count) if canonical_count else 0.0
    ambiguous_content = sum(1 for base in seq if base not in "ACGT") / length if length else 0.0
    return SequenceQuality(
        length=length,
        gc_content=gc_content,
        ambiguous_content=ambiguous_content,
    )


def tree_shrink_qc(
    outdir: Path,
    marker: Marker,
    report: list[BuildReportEntry],
    tree_shrink_config: TreeShrinkConfig,
) -> list[BuildReportEntry]:
    fasta_path = outdir / f"{marker.value}.fasta"
    updated_report = report

    if any(entry.included for entry in report) and fasta_path.exists():
        workdir = outdir / "treeshrink_qc"
        input_fasta = workdir / "input.fasta"
        alignment_path = workdir / "mafft.fasta"
        tree_path = workdir / "iqtree.tree"
        tree_shrink_output_dir = workdir / "treeshrink"
        input_fasta.parent.mkdir(parents=True, exist_ok=True)
        input_fasta.write_text(fasta_path.read_text(encoding="utf-8"), encoding="utf-8")
        AlignmentRunner().align(
            input_fasta,
            alignment_path,
            threads=1,
        )
        TreeRunner().build_tree(
            alignment_path,
            tree_path,
            bootstrap=tree_shrink_config.bootstrap,
        )
        tree_shrink_result = TreeShrinkRunner().detect_outliers(
            tree_path,
            tree_shrink_output_dir,
            quantile=tree_shrink_config.quantile,
            max_removed=tree_shrink_config.max_removed,
        )
        records = [
            record
            for record in SeqIO.parse(str(input_fasta), "fasta")
            if record.id not in tree_shrink_result.removed_taxa
        ]
        with fasta_path.open("w", encoding="utf-8") as handle:
            SeqIO.write(records, handle, "fasta")

        updated_report = []
        for entry in report:
            if entry.included and entry.output_id in tree_shrink_result.removed_taxa:
                metadata = dict(entry.metadata)
                metadata.update(
                    {
                        "tree_shrink_alignment": str(alignment_path),
                        "tree_shrink_tree": str(tree_path),
                        "tree_shrink_output_dir": str(tree_shrink_result.output_dir),
                        "tree_shrink_removed_taxa": str(tree_shrink_result.removed_taxa_path),
                    }
                )
                updated_report.append(
                    replace(
                        entry,
                        included=False,
                        reason="TreeShrink long-branch outlier",
                        output_id=None,
                        metadata=metadata,
                    )
                )
            else:
                updated_report.append(entry)

    (outdir / "build_report.json").write_text(
        json.dumps([asdict(entry) for entry in updated_report], ensure_ascii=False, indent=2)
        + "\n",
        encoding="utf-8",
    )
    return updated_report
