from __future__ import annotations

import io
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from barcode_kit.exceptions import GenBankError
from barcode_kit.models import Marker


__all__ = [
    "AccessionVersion",
    "AnnotationMarkerEvidence",
    "ParsedGenBankRecord",
    "accession_from_record",
    "annotation_marker_evidence",
    "detect_markers",
    "extract_its",
    "extract_marker",
    "format_fasta_record",
    "parse_accession_version",
    "parse_genbank_file",
    "parse_genbank_record",
    "parse_genbank_text",
    "read_single_genbank",
    "source_organism",
    "source_taxon_id",
]


ACCESSION_RE = re.compile(r"^(?P<root>[A-Za-z_]+\d+)(?:\.(?P<version>\d+))?$")


@dataclass(frozen=True)
class AccessionVersion:
    root: str
    version: int

    @property
    def value(self) -> str:
        return f"{self.root}.{self.version}"


@dataclass(frozen=True)
class ParsedGenBankRecord:
    record: SeqRecord
    accession: AccessionVersion
    organism: str
    taxon_id: int | None
    marker_flags: dict[Marker, bool]


@dataclass(frozen=True)
class AnnotationMarkerEvidence:
    contains_marker: bool
    extractable_marker: bool
    sequence: Seq | None = None
    annotation_pattern: str | None = None


@dataclass(frozen=True)
class _ItsFeatureEvidence:
    feature: SeqFeature
    components: frozenset[str]
    start: int
    end: int
    strand: int | None


ITS_FEATURE_TYPES = {"misc_feature", "repeat_region", "rRNA", "misc_RNA"}
ITS_PRIMARY_COMPONENTS = frozenset({"its1", "5_8s", "its2"})
ITS_FLANK_COMPONENTS = frozenset({"ssu", "lsu"})
ITS_QUALIFIER_KEYS = ("product", "note")
ITS1_RE = re.compile(
    r"(?:internal\s+transcribed\s+spacer\s+1(?!\d)|(?<![a-z0-9])its[\s._-]*1(?![a-z0-9]))",
    re.IGNORECASE,
)
ITS2_RE = re.compile(
    r"(?:internal\s+transcribed\s+spacer\s+2(?!\d)|(?<![a-z0-9])its[\s._-]*2(?![a-z0-9]))",
    re.IGNORECASE,
)
RRNA_5_8S_RE = re.compile(
    r"(?:5\s*[\._]\s*8\s*s|5\s*8\s*s)(?:\s+(?:rrna|ribosomal\s+rna))?",
    re.IGNORECASE,
)
SSU_RE = re.compile(r"(?:\bssu\b|\b18\s*s\b|small\s+subunit)", re.IGNORECASE)
LSU_RE = re.compile(r"(?:\blsu\b|\b26\s*s\b|\b28\s*s\b|large\s+subunit)", re.IGNORECASE)


def parse_accession_version(value: str) -> AccessionVersion:
    match = ACCESSION_RE.match(value.strip())
    if not match or not match.group("version"):
        raise GenBankError(f"accession must include version: {value}")
    return AccessionVersion(root=match.group("root"), version=int(match.group("version")))


def read_single_genbank(path: Path) -> SeqRecord:
    records = list(SeqIO.parse(str(path), "genbank"))
    if len(records) != 1:
        raise GenBankError(f"{path} must contain exactly one GenBank record")
    return records[0]


def parse_genbank_text(text: str) -> ParsedGenBankRecord:
    try:
        record = SeqIO.read(io.StringIO(text), "genbank")
    except Exception as error:
        raise GenBankError(f"failed to parse GenBank record: {error}") from error
    return parse_genbank_record(record)


def parse_genbank_file(path: Path) -> ParsedGenBankRecord:
    return parse_genbank_record(read_single_genbank(path))


def parse_genbank_record(record: SeqRecord) -> ParsedGenBankRecord:
    accession = accession_from_record(record)
    organism = str(record.annotations.get("organism") or "").strip()
    if not organism:
        organism = source_organism(record) or "unknown"
    return ParsedGenBankRecord(
        record=record,
        accession=accession,
        organism=organism,
        taxon_id=source_taxon_id(record),
        marker_flags=detect_markers(record),
    )


def accession_from_record(record: SeqRecord) -> AccessionVersion:
    root: str | None = None
    version: int | None = None
    accessions = record.annotations.get("accessions") or []
    if accessions:
        root = str(accessions[0])
    sequence_version = record.annotations.get("sequence_version")
    if sequence_version:
        version = int(sequence_version)
    record_id = str(record.id or "")
    if "." in record_id:
        parsed = parse_accession_version(record_id)
        root = root or parsed.root
        version = version or parsed.version
    if root and version:
        return AccessionVersion(root=root, version=int(version))
    raise GenBankError(f"record missing versioned accession: {record_id}")


def source_organism(record: SeqRecord) -> str | None:
    for feature in record.features:
        if feature.type == "source":
            values = feature.qualifiers.get("organism")
            if values:
                return str(values[0])
    return None


def source_taxon_id(record: SeqRecord) -> int | None:
    for feature in record.features:
        if feature.type != "source":
            continue
        for value in feature.qualifiers.get("db_xref", []):
            if str(value).startswith("taxon:"):
                return int(str(value).split(":", 1)[1])
    return None


def detect_markers(record: SeqRecord) -> dict[Marker, bool]:
    return {marker: annotation_marker_evidence(record, marker).contains_marker for marker in Marker}


def extract_marker(record: SeqRecord, marker: Marker) -> Seq | None:
    return annotation_marker_evidence(record, marker).sequence


def annotation_marker_evidence(record: SeqRecord, marker: Marker) -> AnnotationMarkerEvidence:
    if marker is Marker.ITS:
        return _its_marker_evidence(record)
    if marker is Marker.ITS2:
        return _its2_marker_evidence(record)
    sequence = _extract_coding_marker(record, marker)
    return AnnotationMarkerEvidence(
        contains_marker=sequence is not None,
        extractable_marker=sequence is not None,
        sequence=sequence,
        annotation_pattern="coding_feature" if sequence is not None else None,
    )


def _extract_coding_marker(record: SeqRecord, marker: Marker) -> Seq | None:
    targets = {
        Marker.RBCL: {"rbcl"},
        Marker.MATK: {"matk"},
    }[marker]
    for feature in record.features:
        if feature.type in {"gene", "CDS", "misc_feature", "rRNA", "misc_RNA"} and _qualifier_contains(
            feature.qualifiers, targets
        ):
            return feature.extract(record.seq)
    return None


def extract_its(record: SeqRecord) -> Seq | None:
    return annotation_marker_evidence(record, Marker.ITS).sequence


def _its_marker_evidence(record: SeqRecord) -> AnnotationMarkerEvidence:
    evidence = _its_feature_evidence(record)
    sequence, pattern = _extract_complete_its(record, evidence)
    return AnnotationMarkerEvidence(
        contains_marker=_contains_complete_its(evidence),
        extractable_marker=sequence is not None,
        sequence=sequence,
        annotation_pattern=pattern or _classify_its_pattern(evidence),
    )


def _its2_marker_evidence(record: SeqRecord) -> AnnotationMarkerEvidence:
    evidence = _its_feature_evidence(record)
    sequence = _extract_isolated_its2(record, evidence)
    return AnnotationMarkerEvidence(
        contains_marker=any("its2" in item.components for item in evidence),
        extractable_marker=sequence is not None,
        sequence=sequence,
        annotation_pattern=_classify_its_pattern(evidence),
    )


def _its_feature_evidence(record: SeqRecord) -> list[_ItsFeatureEvidence]:
    evidence: list[_ItsFeatureEvidence] = []
    for feature in record.features:
        if feature.type not in ITS_FEATURE_TYPES:
            continue
        components = _its_components(feature.qualifiers)
        if not components:
            continue
        evidence.append(
            _ItsFeatureEvidence(
                feature=feature,
                components=frozenset(components),
                start=int(feature.location.start),
                end=int(feature.location.end),
                strand=feature.location.strand,
            )
        )
    return evidence


def _its_components(qualifiers: dict[str, Any]) -> set[str]:
    text = _its_qualifier_text(qualifiers)
    components: set[str] = set()
    if ITS1_RE.search(text):
        components.add("its1")
    if ITS2_RE.search(text):
        components.add("its2")
    if RRNA_5_8S_RE.search(text):
        components.add("5_8s")
    if SSU_RE.search(text):
        components.add("ssu")
    if LSU_RE.search(text):
        components.add("lsu")
    return components


def _its_qualifier_text(qualifiers: dict[str, Any]) -> str:
    chunks: list[str] = []
    for key in ITS_QUALIFIER_KEYS:
        values = qualifiers.get(key, [])
        if isinstance(values, str):
            values = [values]
        chunks.extend(str(value) for value in values)
    return " ".join(chunks)


def _contains_complete_its(evidence: list[_ItsFeatureEvidence]) -> bool:
    observed: set[str] = set()
    for item in evidence:
        observed.update(item.components & ITS_PRIMARY_COMPONENTS)
    return ITS_PRIMARY_COMPONENTS <= observed


def _extract_complete_its(
    record: SeqRecord,
    evidence: list[_ItsFeatureEvidence],
) -> tuple[Seq | None, str | None]:
    single_feature = _complete_single_feature(evidence)
    if single_feature is not None:
        return single_feature.extract(record.seq), "complete_single_feature"

    parts = _complete_separate_features(evidence)
    if parts is not None:
        sequence = _extract_complete_its_span(record, parts)
        if sequence is not None:
            return sequence, "complete_separate_features"

    broad_feature = _complete_broad_feature(evidence)
    if broad_feature is not None:
        return broad_feature.extract(record.seq), "broad_ssu_its_lsu"

    return None, None


def _complete_single_feature(evidence: list[_ItsFeatureEvidence]) -> SeqFeature | None:
    for item in evidence:
        if ITS_PRIMARY_COMPONENTS <= item.components and item.components.isdisjoint(ITS_FLANK_COMPONENTS):
            return item.feature
    return None


def _complete_broad_feature(evidence: list[_ItsFeatureEvidence]) -> SeqFeature | None:
    for item in evidence:
        if ITS_PRIMARY_COMPONENTS <= item.components and item.components & ITS_FLANK_COMPONENTS:
            return item.feature
    return None


def _complete_separate_features(
    evidence: list[_ItsFeatureEvidence],
) -> dict[str, _ItsFeatureEvidence] | None:
    parts: dict[str, _ItsFeatureEvidence] = {}
    for item in evidence:
        primary_components = item.components & ITS_PRIMARY_COMPONENTS
        if not primary_components or not item.components.isdisjoint(ITS_FLANK_COMPONENTS):
            continue
        for component in primary_components:
            parts.setdefault(component, item)
    if set(parts) != set(ITS_PRIMARY_COMPONENTS):
        return None
    return parts


def _extract_complete_its_span(
    record: SeqRecord,
    parts: dict[str, _ItsFeatureEvidence],
) -> Seq | None:
    strands = {item.strand for item in parts.values()}
    if len(strands) != 1:
        return None
    strand = next(iter(strands))
    expected_order = ["its2", "5_8s", "its1"] if strand == -1 else ["its1", "5_8s", "its2"]
    ordered_parts = [parts[name] for name in expected_order]
    for left, right in zip(ordered_parts, ordered_parts[1:]):
        if _same_span(left, right):
            continue
        if left.end > right.start:
            return None
    start = min(item.start for item in ordered_parts)
    end = max(item.end for item in ordered_parts)
    sequence = record.seq[start:end]
    if strand == -1:
        sequence = sequence.reverse_complement()
    return sequence


def _same_span(left: _ItsFeatureEvidence, right: _ItsFeatureEvidence) -> bool:
    return left.start == right.start and left.end == right.end and left.strand == right.strand


def _extract_isolated_its2(record: SeqRecord, evidence: list[_ItsFeatureEvidence]) -> Seq | None:
    for item in evidence:
        primary_components = item.components & ITS_PRIMARY_COMPONENTS
        if primary_components == {"its2"} and item.components.isdisjoint(ITS_FLANK_COMPONENTS):
            return item.feature.extract(record.seq)
    for item in evidence:
        if "its2" in item.components:
            return item.feature.extract(record.seq)
    return None


def _classify_its_pattern(evidence: list[_ItsFeatureEvidence]) -> str | None:
    if not evidence:
        return None
    if _complete_single_feature(evidence) is not None:
        return "complete_single_feature"
    if _complete_separate_features(evidence) is not None:
        return "complete_separate_features"
    if _complete_broad_feature(evidence) is not None:
        return "broad_ssu_its_lsu"

    observed: set[str] = set()
    for item in evidence:
        observed.update(item.components)

    if ITS_PRIMARY_COMPONENTS <= observed and observed & ITS_FLANK_COMPONENTS:
        return "broad_ssu_its_lsu"
    if {"5_8s", "its2", "lsu"} <= observed and "its1" not in observed:
        return "partial_5.8s_its2_lsu"
    if "its2" in observed and observed.isdisjoint({"its1", "5_8s", "ssu", "lsu"}):
        return "its2_only"
    if {"its1", "5_8s"} <= observed and "its2" not in observed:
        return "partial_its1_5.8s"
    if observed == {"its1"}:
        return "its1_only"
    return "partial_its_family"


def _qualifier_contains(qualifiers: dict[str, Any], targets: set[str]) -> bool:
    normalized_targets = {target.lower() for target in targets}
    for values in qualifiers.values():
        for value in values:
            text = str(value).lower()
            if any(target in text for target in normalized_targets):
                return True
    return False


def format_fasta_record(identifier: str, sequence: Seq | str, width: int = 80) -> str:
    seq = str(sequence).upper()
    chunks = [seq[index : index + width] for index in range(0, len(seq), width)]
    return f">{identifier}\n" + "\n".join(chunks) + "\n"
