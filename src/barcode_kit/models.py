from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from enum import StrEnum


class Marker(StrEnum):
    ITS = "its"
    ITS2 = "its2"
    MATK = "matk"
    RBCL = "rbcl"

    @property
    def cache_column(self) -> str:
        return {
            Marker.ITS: "has_its",
            Marker.ITS2: "has_its2",
            Marker.MATK: "has_matk",
            Marker.RBCL: "has_rbcl",
        }[self]

    @property
    def ncbi_search_term(self) -> str:
        return {
            Marker.ITS: '"internal transcribed spacer"',
            Marker.ITS2: '"internal transcribed spacer 2"',
            Marker.MATK: "matk",
            Marker.RBCL: "rbcl",
        }[self]

    @property
    def is_coding(self) -> bool:
        return self in {Marker.MATK, Marker.RBCL}


class ItsExtractionMode(StrEnum):
    ITSXRUST = "itsxrust"
    HMM_BLAST = "hmm-blast"
    ANNOTATION = "annotation"


@dataclass(frozen=True)
class TaxonQuery:
    rank: str
    name: str

    def ncbi_term(self) -> str:
        if self.rank == "taxid":
            return f"txid{self.name}[Organism:exp]"
        return f"{self.name}[Organism]"


@dataclass(frozen=True)
class TaxonomyRecord:
    taxon_id: int
    scientific_name: str
    kingdom: str | None = None
    phylum: str | None = None
    class_name: str | None = None
    order: str | None = None
    family: str | None = None
    genus: str | None = None
    species: str | None = None
    infraspecific_rank: str | None = None
    is_hybrid: bool = False
    is_uncertain: bool = False


@dataclass(frozen=True)
class GenBankCacheRecord:
    accession_root: str
    version: int
    accession_version: str
    taxon_id: int
    has_its: bool = False
    has_matk: bool = False
    has_rbcl: bool = False
    has_its2: bool = False
    updated_at: datetime | None = None


@dataclass(frozen=True)
class SequenceQuality:
    length: int
    gc_content: float
    ambiguous_content: float


@dataclass(frozen=True)
class BuildReportEntry:
    accession_version: str
    scientific_name: str | None
    included: bool
    reason: str | None = None
    quality: SequenceQuality | None = None
    output_id: str | None = None
    metadata: dict[str, str | int | float | bool | None] = field(default_factory=dict)
