from __future__ import annotations

import json
import os
import shutil
import tempfile
import uuid
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from barcode_kit.config import TreeShrinkConfig
from barcode_kit.exceptions import BarcodeKitError
from barcode_kit.models import (
    BuildReportEntry,
    Marker,
    SequenceQuality,
    TaxonExclusion,
)
from barcode_kit.phylogeny import (
    AlignmentRunner,
    TreeRunner,
    TreeShrinkRunner,
    TreeShrinkResult,
)


__all__ = [
    "run_qc",
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


def run_qc(
    dataset: Path,
    *,
    exclude: set[TaxonExclusion],
    min_length: int | None,
    max_ambiguous_content: float | None,
    enable_tree_shrink_qc: bool,
    tree_shrink_config: TreeShrinkConfig,
) -> dict[str, Any]:
    if (
        not exclude
        and min_length is None
        and max_ambiguous_content is None
        and not enable_tree_shrink_qc
    ):
        raise BarcodeKitError("select at least one QC option")

    marker, raw_fasta, build_report_path = _load_dataset_manifest(dataset)
    try:
        build_records = json.loads(build_report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise BarcodeKitError(f"invalid build report: {error}") from error
    if not isinstance(build_records, list) or not all(
        isinstance(item, dict) for item in build_records
    ):
        raise BarcodeKitError("invalid build report: expected a list of records")

    build_by_sequence_id: dict[str, dict[str, Any]] = {}
    for build_record in build_records:
        sequence_id = build_record.get("sequence_id")
        if sequence_id is None:
            continue
        if not isinstance(sequence_id, str) or not sequence_id:
            raise BarcodeKitError("invalid build report sequence_id")
        if sequence_id in build_by_sequence_id:
            raise BarcodeKitError(f"duplicate build report sequence_id: {sequence_id}")
        build_by_sequence_id[sequence_id] = build_record

    raw_records = list(SeqIO.parse(str(raw_fasta), "fasta"))

    qc_records: list[dict[str, Any]] = []
    included_records: list[SeqRecord] = []
    seen_ids: set[str] = set()
    for record in raw_records:
        if record.id in seen_ids:
            raise BarcodeKitError(f"duplicate raw FASTA sequence ID: {record.id}")
        seen_ids.add(record.id)
        build_record = build_by_sequence_id.get(record.id)
        if build_record is None:
            raise BarcodeKitError(
                f"raw FASTA sequence ID not found in build report: {record.id}"
            )

        quality = sequence_quality(record.seq)
        reasons: list[str] = []
        if TaxonExclusion.HYBRID in exclude and build_record.get("is_hybrid") is True:
            reasons.append("hybrid_excluded")
        if (
            TaxonExclusion.UNCERTAIN in exclude
            and build_record.get("is_uncertain") is True
        ):
            reasons.append("uncertain_taxon_excluded")
        rank = build_record.get("infraspecific_rank")
        if isinstance(rank, str) and rank:
            if TaxonExclusion.INFRASPECIFIC in exclude:
                reasons.append("infraspecific_taxon_excluded")
            elif any(item.value == rank for item in exclude):
                reasons.append(f"{rank}_excluded")
        if min_length is not None and quality.length < min_length:
            reasons.append("sequence_too_short")
        if (
            max_ambiguous_content is not None
            and quality.ambiguous_content > max_ambiguous_content
        ):
            reasons.append("ambiguous_content_too_high")
        included = not reasons
        if included:
            included_records.append(record)
        qc_records.append(
            {
                "sequence_id": record.id,
                "included": included,
                "length": quality.length,
                "gc_content": quality.gc_content,
                "ambiguous_content": quality.ambiguous_content,
                "reasons": reasons,
            }
        )

    staged_qc = Path(tempfile.mkdtemp(prefix=".qc-stage-", dir=dataset))
    try:
        output_fasta = staged_qc / f"{marker.value}.fasta"
        _write_fasta(output_fasta, included_records)
        if enable_tree_shrink_qc and included_records:
            tree_shrink_result = _run_tree_shrink(
                output_fasta,
                staged_qc / "treeshrink_qc",
                tree_shrink_config,
            )
            removed_taxa = tree_shrink_result.removed_taxa
            for record in qc_records:
                if record["included"] and record["sequence_id"] in removed_taxa:
                    record["included"] = False
                    record["reasons"].append("tree_shrink_long_branch_outlier")
            included_records = [
                record for record in included_records if record.id not in removed_taxa
            ]
            _write_fasta(output_fasta, included_records)

        report = {
            "checks": {
                "exclude": sorted(item.value for item in exclude),
                "min_length": min_length,
                "max_ambiguous_content": max_ambiguous_content,
                "tree_shrink_qc": enable_tree_shrink_qc,
            },
            "records": qc_records,
        }
        (staged_qc / "qc_report.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        _publish_qc_output(dataset, staged_qc)
        return report
    finally:
        if staged_qc.exists():
            shutil.rmtree(staged_qc)


def _load_dataset_manifest(dataset: Path) -> tuple[Marker, Path, Path]:
    manifest_path = dataset / "dataset.json"
    if not manifest_path.is_file():
        raise BarcodeKitError(f"dataset manifest not found: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise BarcodeKitError(f"invalid dataset.json: {error}") from error
    if not isinstance(manifest, dict):
        raise BarcodeKitError("invalid dataset.json: expected an object")
    if manifest.get("format_version") != 1:
        raise BarcodeKitError(
            f"unsupported dataset format_version: {manifest.get('format_version')}"
        )
    try:
        marker = Marker(str(manifest["marker"]))
        raw_fasta = _dataset_file(dataset, manifest["raw_fasta"], "raw FASTA")
        build_report = _dataset_file(
            dataset,
            manifest["build_report"],
            "build report",
        )
    except KeyError as error:
        raise BarcodeKitError(f"invalid dataset.json: missing {error.args[0]}") from error
    except ValueError as error:
        raise BarcodeKitError(f"invalid dataset marker: {manifest.get('marker')}") from error
    return marker, raw_fasta, build_report


def _dataset_file(dataset: Path, value: object, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise BarcodeKitError(f"invalid dataset.json {label} path")
    relative_path = Path(value)
    if relative_path.is_absolute():
        raise BarcodeKitError(f"invalid dataset.json {label} path")
    dataset_root = dataset.resolve()
    path = (dataset / relative_path).resolve()
    if not path.is_relative_to(dataset_root):
        raise BarcodeKitError(f"invalid dataset.json {label} path")
    if not path.is_file():
        raise BarcodeKitError(f"{label} not found: {path}")
    return path


def _write_fasta(path: Path, records: list[SeqRecord]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        SeqIO.write(records, handle, "fasta")


def _publish_qc_output(dataset: Path, staged_qc: Path) -> None:
    qc_dir = dataset / "qc"
    backup = dataset / f".qc-backup-{uuid.uuid4().hex}"
    had_previous_output = qc_dir.exists()
    if had_previous_output:
        os.replace(qc_dir, backup)
    try:
        os.replace(staged_qc, qc_dir)
    except Exception:
        if had_previous_output and backup.exists():
            os.replace(backup, qc_dir)
        raise
    if backup.exists():
        shutil.rmtree(backup)


def _run_tree_shrink(
    fasta_path: Path,
    workdir: Path,
    tree_shrink_config: TreeShrinkConfig,
) -> TreeShrinkResult:
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
    return TreeShrinkRunner().detect_outliers(
        tree_path,
        tree_shrink_output_dir,
        quantile=tree_shrink_config.quantile,
        max_removed=tree_shrink_config.max_removed,
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
        alignment_path = workdir / "mafft.fasta"
        tree_path = workdir / "iqtree.tree"
        input_fasta = workdir / "input.fasta"
        tree_shrink_result = _run_tree_shrink(
            fasta_path,
            workdir,
            tree_shrink_config,
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
