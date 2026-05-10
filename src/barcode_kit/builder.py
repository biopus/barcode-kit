from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from barcode_kit.blast import (
    BlastQuery,
    BlastRescueResult,
    BlastRunner,
    BlastSeed,
    SubprocessBlastRunner,
)
from barcode_kit.config import AppConfig, ensure_app_dirs
from barcode_kit.exceptions import BuildError
from barcode_kit.itsxrust import (
    ItsxrustExtractionResult,
    ItsxrustInput,
    ItsxrustRunner,
    SubprocessItsxrustRunner,
    default_hmm_path,
)
from barcode_kit.models import BuildReportEntry, ItsExtractionMode, Marker, SequenceQuality, TaxonQuery
from barcode_kit.parser import (
    annotation_marker_evidence,
    extract_marker,
    format_fasta_record,
    read_single_genbank,
)
from barcode_kit.storage import Storage
from barcode_kit.validation import sequence_quality


def build_dataset(
    config: AppConfig,
    storage: Storage,
    query: TaxonQuery,
    marker: Marker,
    outdir: Path,
    *,
    min_length: int | None = None,
    max_n_content: float | None = None,
    exclude_hybrid: bool = False,
    exclude_uncertain: bool = False,
    its_extraction_mode: ItsExtractionMode = ItsExtractionMode.HMM_BLAST,
    itsxrust_runner: ItsxrustRunner | None = None,
    its_hmm_path: Path | None = None,
    blast_runner: BlastRunner | None = None,
) -> list[BuildReportEntry]:
    ensure_app_dirs(config)
    storage.initialize()
    outdir.mkdir(parents=True, exist_ok=True)
    candidates = storage.candidate_records(query, marker)
    if not candidates:
        raise BuildError(f"no cached {marker.value} records found for {query.rank} {query.name}")

    fasta_chunks: list[str] = []
    report: list[BuildReportEntry] = []
    itsxrust_runner = itsxrust_runner or SubprocessItsxrustRunner(config=config.itsxrust)
    its_hmm_path = its_hmm_path or default_hmm_path()
    blast_runner = blast_runner or SubprocessBlastRunner(config.blast_rescue)
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
            max_n_content,
            exclude_hybrid,
            exclude_uncertain,
            itsxrust_runner,
            its_hmm_path,
            blast_runner,
        )

    for cache_record, taxonomy in candidates:
        path = config.genbank_cache_dir / f"{cache_record.accession_version}.gb"
        reason: str | None = None
        quality: SequenceQuality | None = None
        sequence = None
        metadata: dict[str, str | int | float | bool | None] = {
            "taxon_id": taxonomy.taxon_id,
            "marker": marker.value,
            "is_hybrid": taxonomy.is_hybrid,
            "is_uncertain": taxonomy.is_uncertain,
        }

        if exclude_hybrid and taxonomy.is_hybrid:
            reason = "hybrid excluded"
        elif exclude_uncertain and taxonomy.is_uncertain:
            reason = "uncertain taxon excluded"
        elif not path.exists():
            reason = "cached GenBank file missing"
        else:
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
                quality = sequence_quality(sequence, marker)
                reason = _filter_reason(quality, marker, min_length, max_n_content)

        included = reason is None and sequence is not None
        output_id = None
        if included:
            output_id = _fasta_id(cache_record.accession_version, taxonomy.scientific_name)
            fasta_chunks.append(format_fasta_record(output_id, sequence))
        report.append(
            BuildReportEntry(
                accession_version=cache_record.accession_version,
                scientific_name=taxonomy.scientific_name,
                included=included,
                reason=reason,
                quality=quality,
                output_id=output_id,
                metadata=metadata,
            )
        )

    (outdir / f"{marker.value}.fasta").write_text("".join(fasta_chunks), encoding="utf-8")
    (outdir / "build_report.json").write_text(
        json.dumps([asdict(entry) for entry in report], ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


@dataclass(frozen=True)
class _PendingBlastRecord:
    index: int
    accession_version: str
    scientific_name: str
    record: SeqRecord
    metadata: dict[str, str | int | float | bool | None]


def _build_its_hmm_blast_dataset(
    config: AppConfig,
    candidates,
    marker: Marker,
    extraction_mode: ItsExtractionMode,
    outdir: Path,
    min_length: int | None,
    max_n_content: float | None,
    exclude_hybrid: bool,
    exclude_uncertain: bool,
    itsxrust_runner: ItsxrustRunner,
    hmm_path: Path,
    blast_runner: BlastRunner,
) -> list[BuildReportEntry]:
    report_slots: list[BuildReportEntry | None] = []
    fasta_sequences: dict[int, Seq] = {}
    itsxrust_records: list[_PendingBlastRecord] = []
    pending: list[_PendingBlastRecord] = []
    seeds: list[BlastSeed] = []

    for cache_record, taxonomy in candidates:
        index = len(report_slots)
        report_slots.append(None)
        accession_version = cache_record.accession_version
        path = config.genbank_cache_dir / f"{accession_version}.gb"
        metadata: dict[str, str | int | float | bool | None] = {
            "taxon_id": taxonomy.taxon_id,
            "marker": marker.value,
            "is_hybrid": taxonomy.is_hybrid,
            "is_uncertain": taxonomy.is_uncertain,
        }

        reason: str | None = None
        quality: SequenceQuality | None = None
        sequence: Seq | None = None
        if exclude_hybrid and taxonomy.is_hybrid:
            reason = "hybrid excluded"
        elif exclude_uncertain and taxonomy.is_uncertain:
            reason = "uncertain taxon excluded"
        elif not path.exists():
            reason = "cached GenBank file missing"
        else:
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

        if reason is None and sequence is None:
            reason = "marker not extracted"
        if reason is None and sequence is not None:
            quality = sequence_quality(sequence, marker)
            reason = _filter_reason(quality, marker, min_length, max_n_content)
            if reason is None:
                seeds.append(BlastSeed(accession_version=accession_version, sequence=sequence))
        entry = _make_report_entry(
            accession_version,
            taxonomy.scientific_name,
            sequence,
            reason,
            quality,
            metadata,
        )
        report_slots[index] = entry
        if entry.included and sequence is not None:
            fasta_sequences[index] = sequence

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
                quality = sequence_quality(sequence, marker)
                reason = _filter_reason(quality, marker, min_length, max_n_content)
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
            if entry.included and sequence is not None:
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
                quality = sequence_quality(result.sequence, marker)
                reason = _filter_reason(quality, marker, min_length, max_n_content)
            entry = _make_report_entry(
                item.accession_version,
                item.scientific_name,
                result.sequence,
                reason,
                quality,
                metadata,
            )
            report_slots[item.index] = entry
            if entry.included and result.sequence is not None:
                fasta_sequences[item.index] = result.sequence

    report = [entry for entry in report_slots if entry is not None]
    fasta_chunks = [
        format_fasta_record(entry.output_id, fasta_sequences[index])
        for index, entry in enumerate(report_slots)
        if entry is not None and entry.included and entry.output_id is not None
    ]
    (outdir / f"{marker.value}.fasta").write_text("".join(fasta_chunks), encoding="utf-8")
    (outdir / "build_report.json").write_text(
        json.dumps([asdict(entry) for entry in report], ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


def _blast_rescue_pending_records(
    pending: list[_PendingBlastRecord],
    seeds: list[BlastSeed],
    marker: Marker,
    blast_runner: BlastRunner,
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


def _hmm_fallback_reason(marker: Marker, result: ItsxrustExtractionResult) -> str:
    if result.fallback_reason is not None:
        return result.fallback_reason
    return "no_anchor_its2" if marker is Marker.ITS2 else "no_anchor_full"


def _filter_reason(
    quality: SequenceQuality,
    marker: Marker,
    min_length: int | None,
    max_n_content: float | None,
) -> str | None:
    if min_length is not None and quality.length < min_length:
        return "sequence shorter than min_length"
    if max_n_content is not None and quality.n_content > max_n_content:
        return "N content above max_n_content"
    if not marker.is_coding:
        return None
    if quality.has_stop_codon:
        return "coding sequence contains stop codon"
    if quality.has_frameshift:
        return "coding sequence length is not divisible by 3"
    return None


def _fasta_id(accession_version: str, scientific_name: str) -> str:
    return f"{accession_version}|{scientific_name.replace(' ', '_')}"
