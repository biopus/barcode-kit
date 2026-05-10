from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from barcode_kit.exceptions import TaxonomyError
from barcode_kit.models import TaxonomyRecord


UNCERTAIN_RE = re.compile(r"\b(sp\.|cf\.|aff\.|unidentified|unknown|uncultured)\b", re.IGNORECASE)
HYBRID_RE = re.compile(r"(×|\bx\b|\bhybrid\b|notho)", re.IGNORECASE)
LINEAGE_RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")


class ETETaxonomyResolver:
    def __init__(self, ncbi: Any | None = None, dbfile: str | Path | None = None):
        self._ncbi = ncbi
        self.dbfile = Path(dbfile).expanduser() if dbfile is not None else None

    def standardize(self, scientific_name: str, taxon_id_hint: int | None = None) -> TaxonomyRecord:
        taxon_id, standardized_name = self._name_to_taxid(scientific_name, taxon_id_hint)
        lineage_by_rank = self._lineage_by_rank(taxon_id)
        resolved_name = lineage_by_rank.get("species") or standardized_name or scientific_name
        return TaxonomyRecord(
            taxon_id=taxon_id,
            scientific_name=resolved_name,
            kingdom=lineage_by_rank.get("kingdom"),
            phylum=lineage_by_rank.get("phylum"),
            class_name=lineage_by_rank.get("class"),
            order=lineage_by_rank.get("order"),
            family=lineage_by_rank.get("family"),
            genus=lineage_by_rank.get("genus"),
            species=_species_epithet(resolved_name),
            is_hybrid=is_hybrid(scientific_name) or is_hybrid(resolved_name),
            is_uncertain=is_uncertain(scientific_name),
        )

    def _name_to_taxid(self, scientific_name: str, taxon_id_hint: int | None) -> tuple[int, str]:
        if taxon_id_hint is not None:
            return taxon_id_hint, self._translated_name(taxon_id_hint) or scientific_name
        try:
            matches = self.ncbi.get_name_translator([scientific_name])
        except Exception as error:
            raise TaxonomyError(f"ETE could not resolve taxon name: {scientific_name}") from error
        taxids = matches.get(scientific_name) or []
        if not taxids:
            raise TaxonomyError(f"ETE could not resolve taxon name: {scientific_name}")
        taxon_id = int(taxids[0])
        return taxon_id, scientific_name

    def _lineage_by_rank(self, taxon_id: int) -> dict[str, str]:
        try:
            lineage = [int(value) for value in self.ncbi.get_lineage(taxon_id)]
        except Exception as error:
            raise TaxonomyError(f"ETE returned empty lineage for {taxon_id}") from error
        if not lineage:
            raise TaxonomyError(f"ETE returned empty lineage for {taxon_id}")
        try:
            ranks_by_taxid = self.ncbi.get_rank(lineage)
            names_by_taxid = self.ncbi.get_taxid_translator(lineage)
        except Exception as error:
            raise TaxonomyError(f"ETE returned invalid lineage metadata for {taxon_id}") from error
        rank_taxids = {rank: lineage_taxid for lineage_taxid, rank in ranks_by_taxid.items()}
        lineage_by_rank = {
            rank: names_by_taxid[lineage_taxid]
            for rank in LINEAGE_RANKS
            if (lineage_taxid := rank_taxids.get(rank)) is not None and lineage_taxid in names_by_taxid
        }
        if not lineage_by_rank:
            raise TaxonomyError(f"ETE returned invalid lineage metadata for {taxon_id}")
        return lineage_by_rank

    def _translated_name(self, taxon_id: int) -> str | None:
        try:
            return self.ncbi.get_taxid_translator([taxon_id]).get(taxon_id)
        except Exception:
            return None

    @property
    def ncbi(self) -> Any:
        if self._ncbi is None:
            try:
                from ete3 import NCBITaxa
            except ImportError as error:
                raise TaxonomyError("ete3 is not installed. Install the ete3 dependency before sync.") from error
            try:
                self._ncbi = NCBITaxa(dbfile=str(self.dbfile) if self.dbfile is not None else None)
            except Exception as error:
                raise TaxonomyError(f"ETE could not initialize the NCBI taxonomy database: {error}") from error
        return self._ncbi


def is_uncertain(name: str) -> bool:
    return bool(UNCERTAIN_RE.search(name))


def is_hybrid(name: str) -> bool:
    return bool(HYBRID_RE.search(name))


def _species_epithet(scientific_name: str) -> str | None:
    parts = scientific_name.split()
    if len(parts) >= 2:
        return parts[1]
    return None
