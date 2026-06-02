from __future__ import annotations

import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from barcode_kit.blast import (
    BlastQuery,
    BlastRescueResult,
    BlastSeed,
    SubprocessBlastRunner,
)
from barcode_kit.config import AppConfig, TreeShrinkConfig, ensure_app_dirs
from barcode_kit.exceptions import BuildError
from barcode_kit.itsxrust import (
    ItsxrustExtractionResult,
    ItsxrustInput,
    SubprocessItsxrustRunner,
    default_hmm_path,
)
from barcode_kit.models import (
    BuildReportEntry,
    ItsExtractionMode,
    Marker,
    SequenceQuality,
    TaxonQuery,
    TaxonomyRecord,
)
from barcode_kit.parser import (
    annotation_marker_evidence,
    extract_marker,
    format_fasta_record,
    read_single_genbank,
)
from barcode_kit.phylogeny import TreeShrinkQcResult, run_tree_shrink_qc
from barcode_kit.storage import Storage
from barcode_kit.validation import sequence_quality


__all__ = ["build_dataset"]


def build_dataset(
    config: AppConfig,
    storage: Storage,
    query: TaxonQuery,
    marker: Marker,
    outdir: Path,
    *,
    min_length: int | None = None,
    max_ambiguous_content: float | None = None,
    exclude_hybrid: bool = False,
    exclude_uncertain: bool = False,
    its_extraction_mode: ItsExtractionMode = ItsExtractionMode.HMM_BLAST,
    enable_tree_shrink_qc: bool = False,
) -> list[BuildReportEntry]:
    ensure_app_dirs(config)
    storage.initialize()
    outdir.mkdir(parents=True, exist_ok=True)
    candidates = storage.candidate_records(query, marker)
    if not candidates:
        raise BuildError(f"no cached {marker.value} records found for {query.rank} {query.name}")

    report: list[BuildReportEntry] = []
    fasta_sequences: list[Seq | None] = []
    if marker in {Marker.ITS, Marker.ITS2} and its_extraction_mode in {
        ItsExtractionMode.HMM_BLAST,
        ItsExtractionMode.ITSXRUST,
    }:
        return _build_its_hmm_blast_dataset(
            config,
            candidates,
            marker,
            its_extraction_mode,
            outdir,
            min_length,
            max_ambiguous_content,
            exclude_hybrid,
            exclude_uncertain,
            SubprocessItsxrustRunner(config=config.itsxrust),
            default_hmm_path(),
            SubprocessBlastRunner(config.blast_rescue),
            enable_tree_shrink_qc,
        )

    for index, (cache_record, taxonomy) in enumerate(candidates):
        path = config.genbank_cache_dir / f"{cache_record.accession_version}.gb"
        reason: str | None = None
        quality: SequenceQuality | None = None
        sequence: Seq | None = None
        metadata = _candidate_metadata(marker, taxonomy)
        reason = _excluded_candidate_reason(
            path,
            taxonomy,
            exclude_hybrid=exclude_hybrid,
            exclude_uncertain=exclude_uncertain,
        )

        if reason is None:
            try:
                record = read_single_genbank(path)
                if marker in {Marker.ITS, Marker.ITS2}:
                    annotation = annotation_marker_evidence(record, marker)
                    sequence = annotation.sequence
                    metadata.update(
                        {
                            "extraction_mode": ItsExtractionMode.ANNOTATION.value,
                            "annotation_pattern": annotation.annotation_pattern,
                            "annotation_contains_marker": annotation.contains_marker,
                            "annotation_extractable_marker": annotation.extractable_marker,
                            "extraction_backend": "annotation",
                            "fallback_reason": None,
                        }
                    )
                else:
                    sequence = extract_marker(record, marker)
            except Exception as error:
                reason = f"GenBank parse failed: {error}"
            if reason is None and sequence is None:
                reason = "marker not extracted"
            if reason is None and sequence is not None:
                quality, reason = _quality_filter_reason(
                    sequence,
                    min_length,
                    max_ambiguous_content,
                )

        entry = _make_report_entry(
            cache_record.accession_version,
            taxonomy.scientific_name,
            sequence,
            reason,
            quality,
            metadata,
        )
        report.append(entry)
        fasta_sequences.append(sequence if entry.included else None)

    return _finalize_build_outputs(
        outdir,
        marker,
        report,
        fasta_sequences,
        config.tree_shrink_qc,
        enable_tree_shrink_qc,
    )


@dataclass(frozen=True)
class _PendingBlastRecord:
    index: int
    accession_version: str
    scientific_name: str
    record: SeqRecord
    metadata: dict[str, str | int | float | bool | None]


def _candidate_metadata(
    marker: Marker,
    taxonomy: TaxonomyRecord,
) -> dict[str, str | int | float | bool | None]:
    return {
        "taxon_id": taxonomy.taxon_id,
        "marker": marker.value,
        "is_hybrid": taxonomy.is_hybrid,
        "is_uncertain": taxonomy.is_uncertain,
    }


def _excluded_candidate_reason(
    path: Path,
    taxonomy: TaxonomyRecord,
    *,
    exclude_hybrid: bool,
    exclude_uncertain: bool,
) -> str | None:
    if exclude_hybrid and taxonomy.is_hybrid:
        return "hybrid excluded"
    if exclude_uncertain and taxonomy.is_uncertain:
        return "uncertain taxon excluded"
    if not path.exists():
        return "cached GenBank file missing"
    return None


def _build_its_hmm_blast_dataset(
    config: AppConfig,
    candidates,
    marker: Marker,
    extraction_mode: ItsExtractionMode,
    outdir: Path,
    min_length: int | None,
    max_ambiguous_content: float | None,
    exclude_hybrid: bool,
    exclude_uncertain: bool,
    itsxrust_runner: SubprocessItsxrustRunner,
    hmm_path: Path,
    blast_runner: SubprocessBlastRunner,
    enable_tree_shrink_qc: bool,
) -> list[BuildReportEntry]:
    report_slots: list[BuildReportEntry | None] = []
    fasta_sequences: list[Seq | None] = []
    itsxrust_records: list[_PendingBlastRecord] = []
    pending: list[_PendingBlastRecord] = []
    seeds: list[BlastSeed] = []

    for cache_record, taxonomy in candidates:
        index = len(report_slots)
        report_slots.append(None)
        fasta_sequences.append(None)
        accession_version = cache_record.accession_version
        path = config.genbank_cache_dir / f"{accession_version}.gb"
        metadata = _candidate_metadata(marker, taxonomy)

        reason = _excluded_candidate_reason(
            path,
            taxonomy,
            exclude_hybrid=exclude_hybrid,
            exclude_uncertain=exclude_uncertain,
        )
        if reason is None:
            try:
                record = read_single_genbank(path)
                annotation = annotation_marker_evidence(record, marker)
                metadata.update(
                    {
                        "extraction_mode": extraction_mode.value,
                        "annotation_pattern": annotation.annotation_pattern,
                        "annotation_contains_marker": annotation.contains_marker,
                        "annotation_extractable_marker": annotation.extractable_marker,
                        "extraction_backend": None,
                        "fallback_reason": None,
                    }
                )
                itsxrust_records.append(
                    _PendingBlastRecord(
                        index=index,
                        accession_version=accession_version,
                        scientific_name=taxonomy.scientific_name,
                        record=record,
                        metadata=metadata,
                    )
                )
                continue
            except Exception as error:
                reason = f"GenBank parse failed: {error}"

        report_slots[index] = _make_report_entry(
            accession_version,
            taxonomy.scientific_name,
            None,
            reason,
            None,
            metadata,
        )

    if itsxrust_records:
        hmm_results = itsxrust_runner.extract_many(
            [
                ItsxrustInput(
                    accession_version=item.accession_version,
                    record=item.record,
                )
                for item in itsxrust_records
            ],
            marker,
            hmm_path,
        )
        for item in itsxrust_records:
            hmm_result = hmm_results.get(
                item.accession_version,
                ItsxrustExtractionResult(
                    sequence=None,
                    fallback_reason="no_anchor_its2" if marker is Marker.ITS2 else "no_anchor_full",
                ),
            )
            metadata = dict(item.metadata)
            sequence = hmm_result.sequence
            reason = None
            quality = None
            if sequence is not None:
                metadata["extraction_backend"] = "itsxrust"
                metadata["fallback_reason"] = None
                quality, reason = _quality_filter_reason(
                    sequence,
                    min_length,
                    max_ambiguous_content,
                )
                if reason is None:
                    seeds.append(BlastSeed(accession_version=item.accession_version, sequence=sequence))
            else:
                metadata["hmm_fallback_reason"] = _hmm_fallback_reason(marker, hmm_result)
                metadata["fallback_reason"] = metadata["hmm_fallback_reason"]
                if extraction_mode is ItsExtractionMode.HMM_BLAST:
                    pending.append(
                        _PendingBlastRecord(
                            index=item.index,
                            accession_version=item.accession_version,
                            scientific_name=item.scientific_name,
                            record=item.record,
                            metadata=metadata,
                        )
                    )
                    continue
                metadata["extraction_backend"] = "itsxrust"
                reason = "marker not extracted"
            entry = _make_report_entry(
                item.accession_version,
                item.scientific_name,
                sequence,
                reason,
                quality,
                metadata,
            )
            report_slots[item.index] = entry
            if entry.included:
                fasta_sequences[item.index] = sequence

    if pending:
        blast_results = _blast_rescue_pending_records(pending, seeds, marker, blast_runner)
        for item in pending:
            result = blast_results.get(
                item.accession_version,
                BlastRescueResult(sequence=None, fallback_reason="no_blast_hit"),
            )
            metadata = dict(item.metadata)
            metadata["extraction_backend"] = "blastn"
            metadata["fallback_reason"] = result.fallback_reason
            metadata.update(result.metadata)
            reason = None
            quality = None
            if result.sequence is None:
                reason = "marker not extracted"
            else:
                quality, reason = _quality_filter_reason(
                    result.sequence,
                    min_length,
                    max_ambiguous_content,
                )
            entry = _make_report_entry(
                item.accession_version,
                item.scientific_name,
                result.sequence,
                reason,
                quality,
                metadata,
            )
            report_slots[item.index] = entry
            if entry.included:
                fasta_sequences[item.index] = result.sequence

    report: list[BuildReportEntry] = []
    included_sequences: list[Seq | None] = []
    for index, entry in enumerate(report_slots):
        if entry is None:
            continue
        report.append(entry)
        included_sequences.append(fasta_sequences[index])
    return _finalize_build_outputs(
        outdir,
        marker,
        report,
        included_sequences,
        config.tree_shrink_qc,
        enable_tree_shrink_qc,
    )


def _blast_rescue_pending_records(
    pending: list[_PendingBlastRecord],
    seeds: list[BlastSeed],
    marker: Marker,
    blast_runner: SubprocessBlastRunner,
) -> dict[str, BlastRescueResult]:
    if not seeds:
        return {
            item.accession_version: BlastRescueResult(
                sequence=None,
                fallback_reason="no_blast_seed",
            )
            for item in pending
        }
    return blast_runner.rescue(
        [
            BlastQuery(
                accession_version=item.accession_version,
                sequence=item.record.seq,
            )
            for item in pending
        ],
        seeds,
        marker,
    )


def _make_report_entry(
    accession_version: str,
    scientific_name: str,
    sequence: Seq | None,
    reason: str | None,
    quality: SequenceQuality | None,
    metadata: dict[str, str | int | float | bool | None],
) -> BuildReportEntry:
    included = reason is None and sequence is not None
    output_id = _fasta_id(accession_version, scientific_name) if included else None
    return BuildReportEntry(
        accession_version=accession_version,
        scientific_name=scientific_name,
        included=included,
        reason=reason,
        quality=quality,
        output_id=output_id,
        metadata=metadata,
    )


def _quality_filter_reason(
    sequence: Seq,
    min_length: int | None,
    max_ambiguous_content: float | None,
) -> tuple[SequenceQuality, str | None]:
    quality = sequence_quality(sequence)
    return quality, _filter_reason(quality, min_length, max_ambiguous_content)


def _write_build_outputs(
    outdir: Path,
    marker: Marker,
    report: list[BuildReportEntry],
    fasta_chunks: list[str],
) -> None:
    (outdir / f"{marker.value}.fasta").write_text("".join(fasta_chunks), encoding="utf-8")
    (outdir / "build_report.json").write_text(
        json.dumps([asdict(entry) for entry in report], ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def _finalize_build_outputs(
    outdir: Path,
    marker: Marker,
    report: list[BuildReportEntry],
    fasta_sequences: list[Seq | None],
    tree_shrink_config: TreeShrinkConfig,
    enable_tree_shrink_qc: bool,
) -> list[BuildReportEntry]:
    fasta_chunks = [
        format_fasta_record(entry.output_id, sequence)
        for entry, sequence in zip(report, fasta_sequences)
        if entry.included and entry.output_id is not None and sequence is not None
    ]
    if not enable_tree_shrink_qc or not fasta_chunks:
        _write_build_outputs(outdir, marker, report, fasta_chunks)
        return report

    workdir = outdir / "treeshrink_qc"
    input_fasta = workdir / "input.fasta"
    input_fasta.parent.mkdir(parents=True, exist_ok=True)
    input_fasta.write_text("".join(fasta_chunks), encoding="utf-8")
    qc_result = run_tree_shrink_qc(
        input_fasta,
        outdir / f"{marker.value}.fasta",
        workdir,
        quantile=tree_shrink_config.quantile,
        bootstrap=tree_shrink_config.bootstrap,
        max_removed=tree_shrink_config.max_removed,
    )
    updated_report = _mark_tree_shrink_outliers(report, qc_result)
    _write_build_report(outdir, updated_report)
    return updated_report


def _mark_tree_shrink_outliers(
    report: list[BuildReportEntry],
    qc_result: TreeShrinkQcResult,
) -> list[BuildReportEntry]:
    updated_report: list[BuildReportEntry] = []
    for entry in report:
        if entry.included and entry.output_id in qc_result.removed_taxa:
            metadata = dict(entry.metadata)
            metadata.update(
                {
                    "tree_shrink_alignment": str(qc_result.alignment_path),
                    "tree_shrink_tree": str(qc_result.tree_path),
                    "tree_shrink_output_dir": str(qc_result.tree_shrink_output_dir),
                    "tree_shrink_removed_taxa": str(qc_result.removed_taxa_path),
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
    return updated_report


def _write_build_report(outdir: Path, report: list[BuildReportEntry]) -> None:
    (outdir / "build_report.json").write_text(
        json.dumps([asdict(entry) for entry in report], ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def _hmm_fallback_reason(marker: Marker, result: ItsxrustExtractionResult) -> str:
    if result.fallback_reason is not None:
        return result.fallback_reason
    return "no_anchor_its2" if marker is Marker.ITS2 else "no_anchor_full"


def _filter_reason(
    quality: SequenceQuality,
    min_length: int | None,
    max_ambiguous_content: float | None,
) -> str | None:
    if min_length is not None and quality.length < min_length:
        return "sequence shorter than min_length"
    if max_ambiguous_content is not None and quality.ambiguous_content > max_ambiguous_content:
        return "ambiguous base content above max_ambiguous_content"
    return None


def _fasta_id(accession_version: str, scientific_name: str) -> str:
    return f"{accession_version}|{scientific_name.replace(' ', '_')}"
