from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from Bio.Seq import Seq

from barcode_kit.config import BlastRescueConfig, BlastRescueMarkerConfig
from barcode_kit.models import Marker
from barcode_kit.parser import format_fasta_record


__all__ = [
    "BlastHit",
    "BlastQuery",
    "BlastRescueDecision",
    "BlastRescueResult",
    "BlastSeed",
    "SubprocessBlastRunner",
    "extract_query_span",
    "parse_blast_tabular",
    "select_blast_rescue_hit",
]


MAKEBLASTDB_COMMAND = "makeblastdb"
BLASTN_COMMAND = "blastn"
BLASTN_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"


@dataclass(frozen=True)
class BlastSeed:
    accession_version: str
    sequence: Seq | str


@dataclass(frozen=True)
class BlastQuery:
    accession_version: str
    sequence: Seq | str


@dataclass(frozen=True)
class BlastHit:
    query_id: str
    subject_id: str
    identity: float
    alignment_length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    mismatches: int = 0
    gap_opens: int = 0


@dataclass(frozen=True)
class BlastRescueDecision:
    accepted: bool
    hit: BlastHit | None = None
    fallback_reason: str | None = None
    subject_coverage: float | None = None
    query_start: int | None = None
    query_end: int | None = None
    strand: str | None = None


@dataclass(frozen=True)
class BlastRescueResult:
    sequence: Seq | None
    fallback_reason: str | None
    metadata: dict[str, str | int | float | bool | None] = field(default_factory=dict)


def select_blast_rescue_hit(
    hits: Iterable[BlastHit],
    marker: Marker,
    seed_lengths: dict[str, int],
    config: BlastRescueConfig | None = None,
) -> BlastRescueDecision:
    config = config or BlastRescueConfig()
    hit_list = list(hits)
    if not hit_list:
        return BlastRescueDecision(accepted=False, fallback_reason="no_blast_hit")

    evaluated = [_evaluate_hit(hit, marker, seed_lengths, config) for hit in hit_list]
    accepted = [decision for decision in evaluated if decision.accepted]
    if not accepted:
        return max(
            evaluated,
            key=lambda decision: decision.hit.bitscore if decision.hit is not None else 0.0,
        )

    accepted.sort(
        key=lambda decision: (
            decision.hit.bitscore if decision.hit is not None else 0.0,
            decision.hit.identity if decision.hit is not None else 0.0,
            decision.subject_coverage or 0.0,
        ),
        reverse=True,
    )
    best = accepted[0]
    for other in accepted[1:]:
        if best.hit is None or other.hit is None:
            continue
        if other.hit.bitscore < best.hit.bitscore * config.ambiguous_bitscore_ratio:
            continue
        if _query_overlap_ratio(best.hit, other.hit) < config.ambiguous_overlap_ratio:
            return BlastRescueDecision(
                accepted=False,
                hit=best.hit,
                fallback_reason="ambiguous_blast_hit",
                subject_coverage=best.subject_coverage,
                query_start=best.query_start,
                query_end=best.query_end,
                strand=best.strand,
            )
    return best


def extract_query_span(sequence: Seq | str, hit: BlastHit) -> Seq:
    start = min(hit.qstart, hit.qend) - 1
    end = max(hit.qstart, hit.qend)
    extracted = Seq(str(sequence)[start:end])
    if _hit_strand(hit) == "-":
        return extracted.reverse_complement()
    return extracted


def parse_blast_tabular(text: str) -> list[BlastHit]:
    hits: list[BlastHit] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 12:
            continue
        hits.append(
            BlastHit(
                query_id=fields[0],
                subject_id=fields[1],
                identity=float(fields[2]),
                alignment_length=int(fields[3]),
                mismatches=int(fields[4]),
                gap_opens=int(fields[5]),
                qstart=int(fields[6]),
                qend=int(fields[7]),
                sstart=int(fields[8]),
                send=int(fields[9]),
                evalue=float(fields[10]),
                bitscore=float(fields[11]),
            )
        )
    return hits


class SubprocessBlastRunner:
    def __init__(self, config: BlastRescueConfig | None = None):
        self.config = config or BlastRescueConfig()

    def rescue(
        self,
        failed_records: Sequence[BlastQuery],
        seeds: Sequence[BlastSeed],
        marker: Marker,
    ) -> dict[str, BlastRescueResult]:
        if not failed_records:
            return {}
        if not seeds:
            return {
                record.accession_version: BlastRescueResult(
                    sequence=None,
                    fallback_reason="no_blast_seed",
                )
                for record in failed_records
            }
        if shutil.which(MAKEBLASTDB_COMMAND) is None or shutil.which(BLASTN_COMMAND) is None:
            return _failed_results(failed_records, "blast_tool_unavailable")

        with tempfile.TemporaryDirectory(prefix="barcode-kit-blast-") as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            seed_path = temp_dir / "seeds.fasta"
            query_path = temp_dir / "queries.fasta"
            db_path = temp_dir / "seeds_db"
            seed_path.write_text(
                "".join(format_fasta_record(seed.accession_version, seed.sequence) for seed in seeds),
                encoding="utf-8",
            )
            query_path.write_text(
                "".join(format_fasta_record(record.accession_version, record.sequence) for record in failed_records),
                encoding="utf-8",
            )

            make_db = subprocess.run(
                [
                    MAKEBLASTDB_COMMAND,
                    "-in",
                    str(seed_path),
                    "-dbtype",
                    "nucl",
                    "-out",
                    str(db_path),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            if make_db.returncode != 0:
                return _failed_results(failed_records, "makeblastdb_failed")

            blast = subprocess.run(
                [
                    BLASTN_COMMAND,
                    "-query",
                    str(query_path),
                    "-db",
                    str(db_path),
                    "-outfmt",
                    BLASTN_OUTFMT,
                    "-dust",
                    self.config.blastn_dust,
                    "-word_size",
                    str(self.config.word_size),
                    "-evalue",
                    str(self.config.evalue),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            if blast.returncode != 0:
                return _failed_results(failed_records, "blastn_failed")

        grouped: dict[str, list[BlastHit]] = defaultdict(list)
        for hit in parse_blast_tabular(blast.stdout):
            grouped[hit.query_id].append(hit)

        seed_lengths = {seed.accession_version: len(str(seed.sequence)) for seed in seeds}
        results: dict[str, BlastRescueResult] = {}
        for record in failed_records:
            decision = select_blast_rescue_hit(
                grouped.get(record.accession_version, []),
                marker,
                seed_lengths,
                self.config,
            )
            if not decision.accepted or decision.hit is None:
                results[record.accession_version] = BlastRescueResult(
                    sequence=None,
                    fallback_reason=decision.fallback_reason,
                    metadata=_decision_metadata(decision),
                )
                continue
            results[record.accession_version] = BlastRescueResult(
                sequence=extract_query_span(record.sequence, decision.hit),
                fallback_reason=None,
                metadata=_decision_metadata(decision),
            )
        return results


def _evaluate_hit(
    hit: BlastHit,
    marker: Marker,
    seed_lengths: dict[str, int],
    config: BlastRescueConfig,
) -> BlastRescueDecision:
    seed_length = seed_lengths.get(hit.subject_id)
    if seed_length is None or seed_length <= 0:
        return BlastRescueDecision(accepted=False, hit=hit, fallback_reason="blast_seed_unknown")

    subject_start = min(hit.sstart, hit.send)
    subject_end = max(hit.sstart, hit.send)
    subject_coverage = (subject_end - subject_start + 1) / seed_length
    query_span_length = abs(hit.qend - hit.qstart) + 1
    length_ratio = query_span_length / seed_length
    thresholds = _thresholds(marker, config)
    query_start = min(hit.qstart, hit.qend)
    query_end = max(hit.qstart, hit.qend)

    if subject_coverage < thresholds.min_subject_coverage:
        return BlastRescueDecision(
            accepted=False,
            hit=hit,
            fallback_reason="blast_subject_coverage_low",
            subject_coverage=subject_coverage,
            query_start=query_start,
            query_end=query_end,
            strand=_hit_strand(hit),
        )
    if hit.identity / 100 < thresholds.min_identity:
        return BlastRescueDecision(
            accepted=False,
            hit=hit,
            fallback_reason="blast_identity_low",
            subject_coverage=subject_coverage,
            query_start=query_start,
            query_end=query_end,
            strand=_hit_strand(hit),
        )
    endpoint_margin = max(
        config.endpoint_margin_bases,
        seed_length * config.endpoint_margin_fraction,
    )
    if subject_start > endpoint_margin or subject_end < seed_length - endpoint_margin:
        return BlastRescueDecision(
            accepted=False,
            hit=hit,
            fallback_reason="blast_seed_endpoint_miss",
            subject_coverage=subject_coverage,
            query_start=query_start,
            query_end=query_end,
            strand=_hit_strand(hit),
        )
    if (
        length_ratio < thresholds.min_query_length_ratio
        or length_ratio > thresholds.max_query_length_ratio
    ):
        return BlastRescueDecision(
            accepted=False,
            hit=hit,
            fallback_reason="blast_query_length_ratio",
            subject_coverage=subject_coverage,
            query_start=query_start,
            query_end=query_end,
            strand=_hit_strand(hit),
        )
    return BlastRescueDecision(
        accepted=True,
        hit=hit,
        subject_coverage=subject_coverage,
        query_start=query_start,
        query_end=query_end,
        strand=_hit_strand(hit),
    )


def _thresholds(marker: Marker, config: BlastRescueConfig) -> BlastRescueMarkerConfig:
    if marker is Marker.ITS2:
        return config.its2
    return config.its


def _hit_strand(hit: BlastHit) -> str:
    query_forward = hit.qend >= hit.qstart
    subject_forward = hit.send >= hit.sstart
    return "+" if query_forward == subject_forward else "-"


def _query_overlap_ratio(left: BlastHit, right: BlastHit) -> float:
    left_start, left_end = sorted((left.qstart, left.qend))
    right_start, right_end = sorted((right.qstart, right.qend))
    overlap = max(0, min(left_end, right_end) - max(left_start, right_start) + 1)
    shorter = min(left_end - left_start + 1, right_end - right_start + 1)
    if shorter <= 0:
        return 0.0
    return overlap / shorter


def _decision_metadata(decision: BlastRescueDecision) -> dict[str, str | int | float | bool | None]:
    if decision.hit is None:
        return {}
    return {
        "blast_seed_accession": decision.hit.subject_id,
        "blast_identity": decision.hit.identity,
        "blast_subject_coverage": decision.subject_coverage,
        "blast_query_start": decision.query_start,
        "blast_query_end": decision.query_end,
        "blast_bitscore": decision.hit.bitscore,
        "blast_strand": decision.strand,
    }


def _failed_results(
    failed_records: Sequence[BlastQuery],
    fallback_reason: str,
) -> dict[str, BlastRescueResult]:
    return {
        record.accession_version: BlastRescueResult(sequence=None, fallback_reason=fallback_reason)
        for record in failed_records
    }
