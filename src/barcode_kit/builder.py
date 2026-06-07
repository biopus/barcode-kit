from __future__ import annotations

import json
from dataclasses import dataclass
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


__all__ = ["build_dataset"]


def build_dataset(
    config: AppConfig,
    storage: Storage,
    query: TaxonQuery,
    marker: Marker,
    outdir: Path,
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
        )

    if marker.is_coding:
        return _build_annotation_dataset(
            config,
            candidates,
            marker,
            outdir,
        )

    raise BuildError(f"unsupported marker: {marker.value}")


_CandidateRecords = list[tuple[GenBankCacheRecord, TaxonomyRecord]]


def _build_its_dataset(
    config: AppConfig,
    candidates: _CandidateRecords,
    marker: Marker,
    outdir: Path,
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
            "infraspecific_rank": taxonomy.infraspecific_rank,
        }

        reason = _build_error(path)
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
            except Exception:
                reason = "genbank_parse_failed"

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
            if sequence is not None:
                metadata["extraction_backend"] = "itsxrust"
                metadata["fallback_reason"] = None
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
                quality=None,
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
            if sequence is None:
                reason = "marker_not_extracted"
            included = reason is None and sequence is not None
            entry = BuildReportEntry(
                accession_version=item.accession_version,
                scientific_name=item.scientific_name,
                included=included,
                reason=reason,
                quality=None,
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
            if result.sequence is None:
                reason = "marker_not_extracted"
            included = reason is None and result.sequence is not None
            entry = BuildReportEntry(
                accession_version=item.accession_version,
                scientific_name=item.scientific_name,
                included=included,
                reason=reason,
                quality=None,
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
) -> list[BuildReportEntry]:
    report: list[BuildReportEntry] = []
    fasta_sequences: list[Seq | None] = []
    for cache_record, taxonomy in candidates:
        path = config.genbank_cache_dir / f"{cache_record.accession_version}.gb"
        reason: str | None = None
        sequence: Seq | None = None
        metadata = {
            "taxon_id": taxonomy.taxon_id,
            "marker": marker.value,
            "is_hybrid": taxonomy.is_hybrid,
            "is_uncertain": taxonomy.is_uncertain,
            "infraspecific_rank": taxonomy.infraspecific_rank,
            "extraction_backend": "annotation",
        }
        reason = _build_error(path)

        if reason is None:
            try:
                record = read_single_genbank(path)
                sequence = extract_marker(record, marker)
            except Exception:
                reason = "genbank_parse_failed"
            if reason is None and sequence is None:
                reason = "marker_not_extracted"

        included = reason is None and sequence is not None
        entry = BuildReportEntry(
            accession_version=cache_record.accession_version,
            scientific_name=taxonomy.scientific_name,
            included=included,
            reason=reason,
            quality=None,
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


def _build_error(path: Path) -> str | None:
    if not path.exists():
        return "cached_file_missing"
    return None


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
    raw_dir = outdir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    raw_fasta = raw_dir / f"{marker.value}.fasta"
    raw_fasta.write_text("".join(fasta_chunks), encoding="utf-8")
    (outdir / "build_report.json").write_text(
        json.dumps([_build_report_record(entry) for entry in report], ensure_ascii=False, indent=2)
        + "\n",
        encoding="utf-8",
    )
    (outdir / "dataset.json").write_text(
        json.dumps(
            {
                "format_version": 1,
                "marker": marker.value,
                "raw_fasta": str(raw_fasta.relative_to(outdir)),
                "build_report": "build_report.json",
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return report


def _build_report_record(
    entry: BuildReportEntry,
) -> dict[str, str | bool | None]:
    return {
        "accession": entry.accession_version,
        "sequence_id": entry.output_id,
        "scientific_name": entry.scientific_name,
        "infraspecific_rank": _string_metadata(entry, "infraspecific_rank"),
        "is_hybrid": bool(entry.metadata.get("is_hybrid", False)),
        "is_uncertain": bool(entry.metadata.get("is_uncertain", False)),
        "extraction_backend": (
            _string_metadata(entry, "extraction_backend") if entry.included else None
        ),
        "error": entry.reason,
    }


def _string_metadata(entry: BuildReportEntry, key: str) -> str | None:
    value = entry.metadata.get(key)
    return value if isinstance(value, str) else None
