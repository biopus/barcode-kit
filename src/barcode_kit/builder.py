from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from barcode_kit.blast import (
    BlastRunner,
    BlastRecord,
    BlastRescueResult,
)
from barcode_kit.config import AppConfig, ensure_app_dirs
from barcode_kit.exceptions import BuildError
from barcode_kit.itsxrust import (
    ItsxrustRunner,
    ItsxrustExtractionResult,
    ItsxrustInput,
    default_hmm_path,
)
from barcode_kit.models import (
    BuildReportEntry,
    GenBankCacheRecord,
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
) -> list[BuildReportEntry]:
    ensure_app_dirs(config)
    storage.initialize()
    outdir.mkdir(parents=True, exist_ok=True)
    candidates = storage.candidate_records(query, marker)
    if not candidates:
        raise BuildError(f"no cached {marker.value} records found for {query.rank} {query.name}")

    if marker in {Marker.ITS, Marker.ITS2}:
        return _build_its_dataset(
            config,
            candidates,
            marker,
            outdir,
            min_length=min_length,
            max_ambiguous_content=max_ambiguous_content,
            exclude_hybrid=exclude_hybrid,
            exclude_uncertain=exclude_uncertain,
        )

    if marker.is_coding:
        return _build_annotation_dataset(
            config,
            candidates,
            marker,
            outdir,
            min_length=min_length,
            max_ambiguous_content=max_ambiguous_content,
            exclude_hybrid=exclude_hybrid,
            exclude_uncertain=exclude_uncertain,
        )

    raise BuildError(f"unsupported marker: {marker.value}")


_CandidateRecords = list[tuple[GenBankCacheRecord, TaxonomyRecord]]


def _build_its_dataset(
    config: AppConfig,
    candidates: _CandidateRecords,
    marker: Marker,
    outdir: Path,
    *,
    min_length: int | None,
    max_ambiguous_content: float | None,
    exclude_hybrid: bool,
    exclude_uncertain: bool,
) -> list[BuildReportEntry]:
    report_slots: list[BuildReportEntry | None] = []
    fasta_sequences: list[Seq | None] = []
    its_records: list[_ItsCandidate] = []
    failed_itsxrust_records: list[_ItsCandidate] = []
    seeds: list[BlastRecord] = []

    for cache_record, taxonomy in candidates:
        index = len(report_slots)
        report_slots.append(None)
        fasta_sequences.append(None)
        accession_version = cache_record.accession_version
        path = config.genbank_cache_dir / f"{accession_version}.gb"
        metadata = {
            "taxon_id": taxonomy.taxon_id,
            "marker": marker.value,
            "is_hybrid": taxonomy.is_hybrid,
            "is_uncertain": taxonomy.is_uncertain,
        }

        reason = _meta_info_qc(
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
                        "annotation_pattern": annotation.annotation_pattern,
                        "annotation_contains_marker": annotation.contains_marker,
                        "annotation_extractable_marker": annotation.extractable_marker,
                        "extraction_backend": None,
                        "fallback_reason": None,
                    }
                )
                its_records.append(
                    _ItsCandidate(
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

        report_slots[index] = BuildReportEntry(
            accession_version=accession_version,
            scientific_name=taxonomy.scientific_name,
            included=False,
            reason=reason,
            quality=None,
            output_id=None,
            metadata=metadata,
        )

    if its_records:
        hmm_results = ItsxrustRunner(config=config.itsxrust).extract_many(
            [
                ItsxrustInput(
                    accession_version=item.accession_version,
                    record=item.record,
                )
                for item in its_records
            ],
            marker,
            default_hmm_path(),
        )
        for item in its_records:
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
                quality, reason = _seq_quality_qc(
                    sequence,
                    min_length,
                    max_ambiguous_content,
                )
                if reason is None:
                    seeds.append(BlastRecord(accession_version=item.accession_version, sequence=sequence))
            else:
                default_fallback_reason = "no_anchor_its2" if marker is Marker.ITS2 else "no_anchor_full"
                metadata["hmm_fallback_reason"] = hmm_result.fallback_reason or default_fallback_reason
                metadata["fallback_reason"] = metadata["hmm_fallback_reason"]
                failed_itsxrust_records.append(
                    _ItsCandidate(
                        index=item.index,
                        accession_version=item.accession_version,
                        scientific_name=item.scientific_name,
                        record=item.record,
                        metadata=metadata,
                    )
                )
                continue
            included = reason is None and sequence is not None
            entry = BuildReportEntry(
                accession_version=item.accession_version,
                scientific_name=item.scientific_name,
                included=included,
                reason=reason,
                quality=quality,
                output_id=f"{item.accession_version}|{item.scientific_name.replace(' ', '_')}"
                if included
                else None,
                metadata=metadata,
            )
            report_slots[item.index] = entry
            if entry.included:
                fasta_sequences[item.index] = sequence

    if failed_itsxrust_records and len(failed_itsxrust_records) == len(its_records):
        for item in failed_itsxrust_records:
            annotation = annotation_marker_evidence(item.record, marker)
            metadata = dict(item.metadata)
            metadata.update(
                {
                    "annotation_pattern": annotation.annotation_pattern,
                    "annotation_contains_marker": annotation.contains_marker,
                    "annotation_extractable_marker": annotation.extractable_marker,
                    "extraction_backend": "annotation",
                    "fallback_reason": None,
                }
            )
            sequence = annotation.sequence
            reason = None
            quality = None
            if sequence is None:
                reason = "marker not extracted"
            else:
                quality, reason = _seq_quality_qc(
                    sequence,
                    min_length,
                    max_ambiguous_content,
                )
            included = reason is None and sequence is not None
            entry = BuildReportEntry(
                accession_version=item.accession_version,
                scientific_name=item.scientific_name,
                included=included,
                reason=reason,
                quality=quality,
                output_id=f"{item.accession_version}|{item.scientific_name.replace(' ', '_')}"
                if included
                else None,
                metadata=metadata,
            )
            report_slots[item.index] = entry
            if entry.included:
                fasta_sequences[item.index] = sequence

    elif failed_itsxrust_records:
        if seeds:
            blast_results = BlastRunner(config.blast_rescue).rescue(
                [
                    BlastRecord(
                        accession_version=item.accession_version,
                        sequence=item.record.seq,
                    )
                    for item in failed_itsxrust_records
                ],
                seeds,
                marker,
            )
        else:
            blast_results = {
                item.accession_version: BlastRescueResult(
                    sequence=None,
                    fallback_reason="no_blast_seed",
                )
                for item in failed_itsxrust_records
            }
        for item in failed_itsxrust_records:
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
                reason = "low quality sequence"
            else:
                quality, reason = _seq_quality_qc(
                    result.sequence,
                    min_length,
                    max_ambiguous_content,
                )
            included = reason is None and result.sequence is not None
            entry = BuildReportEntry(
                accession_version=item.accession_version,
                scientific_name=item.scientific_name,
                included=included,
                reason=reason,
                quality=quality,
                output_id=f"{item.accession_version}|{item.scientific_name.replace(' ', '_')}"
                if included
                else None,
                metadata=metadata,
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
    )


def _build_annotation_dataset(
    config: AppConfig,
    candidates: _CandidateRecords,
    marker: Marker,
    outdir: Path,
    *,
    min_length: int | None,
    max_ambiguous_content: float | None,
    exclude_hybrid: bool,
    exclude_uncertain: bool,
) -> list[BuildReportEntry]:
    report: list[BuildReportEntry] = []
    fasta_sequences: list[Seq | None] = []
    for cache_record, taxonomy in candidates:
        path = config.genbank_cache_dir / f"{cache_record.accession_version}.gb"
        reason: str | None = None
        quality: SequenceQuality | None = None
        sequence: Seq | None = None
        metadata = {
            "taxon_id": taxonomy.taxon_id,
            "marker": marker.value,
            "is_hybrid": taxonomy.is_hybrid,
            "is_uncertain": taxonomy.is_uncertain,
        }
        reason = _meta_info_qc(
            path,
            taxonomy,
            exclude_hybrid=exclude_hybrid,
            exclude_uncertain=exclude_uncertain,
        )

        if reason is None:
            try:
                record = read_single_genbank(path)
                sequence = extract_marker(record, marker)
            except Exception as error:
                reason = f"GenBank parse failed: {error}"
            if reason is None and sequence is None:
                reason = "marker not extracted"
            if reason is None and sequence is not None:
                quality, reason = _seq_quality_qc(
                    sequence,
                    min_length,
                    max_ambiguous_content,
                )

        included = reason is None and sequence is not None
        entry = BuildReportEntry(
            accession_version=cache_record.accession_version,
            scientific_name=taxonomy.scientific_name,
            included=included,
            reason=reason,
            quality=quality,
            output_id=f"{cache_record.accession_version}|{taxonomy.scientific_name.replace(' ', '_')}"
            if included
            else None,
            metadata=metadata,
        )
        report.append(entry)
        fasta_sequences.append(sequence if entry.included else None)

    return _finalize_build_outputs(
        outdir,
        marker,
        report,
        fasta_sequences,
    )


@dataclass(frozen=True)
class _ItsCandidate:
    index: int
    accession_version: str
    scientific_name: str
    record: SeqRecord
    metadata: dict[str, str | int | float | bool | None]


def _meta_info_qc(
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


def _seq_quality_qc(
    sequence: Seq,
    min_length: int | None,
    max_ambiguous_content: float | None,
) -> tuple[SequenceQuality, str | None]:
    quality = sequence_quality(sequence)
    reason = None
    if min_length is not None and quality.length < min_length:
        reason = "sequence shorter than min_length"
    if max_ambiguous_content is not None and quality.ambiguous_content > max_ambiguous_content:
        reason = "ambiguous base content above max_ambiguous_content"
    return quality, reason


def _finalize_build_outputs(
    outdir: Path,
    marker: Marker,
    report: list[BuildReportEntry],
    fasta_sequences: list[Seq | None],
) -> list[BuildReportEntry]:
    fasta_chunks = [
        format_fasta_record(entry.output_id, sequence)
        for entry, sequence in zip(report, fasta_sequences)
        if entry.included and entry.output_id is not None and sequence is not None
    ]
    (outdir / f"{marker.value}.fasta").write_text("".join(fasta_chunks), encoding="utf-8")
    (outdir / "build_report.json").write_text(
        json.dumps([asdict(entry) for entry in report], ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report
