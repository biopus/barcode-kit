from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from barcode_kit.config import ItsxrustConfig
from barcode_kit.models import Marker
from barcode_kit.parser import format_fasta_record


__all__ = [
    "ItsxrustExtractionResult",
    "ItsxrustInput",
    "SubprocessItsxrustRunner",
    "default_hmm_path",
]


@dataclass(frozen=True)
class ItsxrustExtractionResult:
    sequence: Seq | None
    fallback_reason: str | None = None


@dataclass(frozen=True)
class ItsxrustInput:
    accession_version: str
    record: SeqRecord


class SubprocessItsxrustRunner:
    def __init__(self, batch_size: int = 1000, config: ItsxrustConfig | None = None):
        self.batch_size = batch_size
        self.config = config or ItsxrustConfig()

    def extract_many(
        self,
        records: Sequence[ItsxrustInput],
        marker: Marker,
        hmm_path: Path,
    ) -> dict[str, ItsxrustExtractionResult]:
        if not records:
            return {}
        if marker not in {Marker.ITS, Marker.ITS2}:
            return _results_for(records, "unsupported_marker")
        if not hmm_path.exists():
            return _results_for(records, "hmm_missing")
        if shutil.which("itsxrust") is None or shutil.which("nhmmer") is None:
            return _results_for(records, "tool_unavailable")

        region = "its2" if marker is Marker.ITS2 else "full"
        no_anchor_reason = "no_anchor_its2" if marker is Marker.ITS2 else "no_anchor_full"
        results: dict[str, ItsxrustExtractionResult] = {}
        for batch in _chunks(records, self.batch_size):
            results.update(self._extract_batch(batch, hmm_path, region, no_anchor_reason))
        return results

    def _extract_batch(
        self,
        records: Sequence[ItsxrustInput],
        hmm_path: Path,
        region: str,
        no_anchor_reason: str,
    ) -> dict[str, ItsxrustExtractionResult]:
        with tempfile.TemporaryDirectory(prefix="barcode-kit-itsxrust-") as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / "input.fasta"
            output_path = temp_dir / "output.fasta"
            input_path.write_text(
                "".join(
                    format_fasta_record(item.accession_version, item.record.seq)
                    for item in records
                ),
                encoding="utf-8",
            )
            command = [
                "itsxrust",
                "extract",
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                "--hmm",
                str(hmm_path),
                "--region",
                region,
                "--output-format",
                "fasta",
                "--hmmer-cpu",
                "1",
                "--inc-e",
                str(self.config.inc_e),
                "--min-anchor-score",
                str(self.config.min_anchor_score),
                "--max-per-anchor",
                str(self.config.max_per_anchor),
                "--max-anchor-evalue",
                str(self.config.max_anchor_evalue),
            ]
            try:
                completed = subprocess.run(
                    command,
                    check=False,
                    capture_output=True,
                    text=True,
                )
            except FileNotFoundError:
                return _results_for(records, "tool_unavailable")
            if completed.returncode != 0:
                return _results_for(records, "hmm_failed")
            if not output_path.exists():
                return _results_for(records, no_anchor_reason)
            output_records = list(SeqIO.parse(str(output_path), "fasta"))
            output_by_id = {record.id.split('|')[0]: record.seq for record in output_records}
            return {
                item.accession_version: ItsxrustExtractionResult(
                    sequence=output_by_id.get(item.accession_version),
                    fallback_reason=None if item.accession_version in output_by_id else no_anchor_reason,
                )
                for item in records
            }


def _chunks(records: Sequence[ItsxrustInput], size: int) -> list[Sequence[ItsxrustInput]]:
    chunk_size = max(1, size)
    return [records[index : index + chunk_size] for index in range(0, len(records), chunk_size)]


def _results_for(
    records: Sequence[ItsxrustInput],
    fallback_reason: str,
) -> dict[str, ItsxrustExtractionResult]:
    return {
        item.accession_version: ItsxrustExtractionResult(sequence=None, fallback_reason=fallback_reason)
        for item in records
    }


def default_hmm_path() -> Path:
    return Path(__file__).resolve().parents[2] / "hmm" / "T.hmm"
