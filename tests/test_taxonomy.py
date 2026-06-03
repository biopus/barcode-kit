from __future__ import annotations

import pytest

from barcode_kit.exceptions import TaxonomyError
from barcode_kit.taxonomy import ETETaxonomyResolver


class FakeNCBITaxa:
    def __init__(self):
        self.name_requests: list[list[str]] = []
        self.lineage_requests: list[int] = []
        self.rank_requests: list[list[int]] = []
        self.taxid_requests: list[list[int]] = []

    def get_name_translator(self, names: list[str]) -> dict[str, list[int]]:
        self.name_requests.append(names)
        return {"Iris japonica": [12345]}

    def get_lineage(self, taxid: int) -> list[int]:
        self.lineage_requests.append(taxid)
        return [33090, 35493, 3398, 73496, 58920, 26379, 12345]

    def get_rank(self, taxids: list[int]) -> dict[int, str]:
        self.rank_requests.append(taxids)
        return {
            33090: "kingdom",
            35493: "phylum",
            3398: "class",
            73496: "order",
            58920: "family",
            26379: "genus",
            12345: "species",
        }

    def get_taxid_translator(self, taxids: list[int]) -> dict[int, str]:
        self.taxid_requests.append(taxids)
        return {
            33090: "Viridiplantae",
            35493: "Streptophyta",
            3398: "Magnoliopsida",
            73496: "Asparagales",
            58920: "Iridaceae",
            26379: "Iris",
            12345: "Iris japonica",
        }


def test_ete_resolver_uses_name_translation_and_lineage_ranks():
    ncbi = FakeNCBITaxa()

    record = ETETaxonomyResolver(ncbi=ncbi).standardize("Iris japonica")

    assert ncbi.name_requests == [["Iris japonica"]]
    assert ncbi.lineage_requests == [12345]
    assert ncbi.rank_requests == [[33090, 35493, 3398, 73496, 58920, 26379, 12345]]
    assert ncbi.taxid_requests == [[33090, 35493, 3398, 73496, 58920, 26379, 12345]]
    assert record.taxon_id == 12345
    assert record.scientific_name == "Iris japonica"
    assert record.kingdom == "Viridiplantae"
    assert record.phylum == "Streptophyta"
    assert record.class_name == "Magnoliopsida"
    assert record.order == "Asparagales"
    assert record.family == "Iridaceae"
    assert record.genus == "Iris"
    assert record.species == "japonica"


def test_ete_resolver_prefers_taxon_id_hint():
    ncbi = FakeNCBITaxa()

    record = ETETaxonomyResolver(ncbi=ncbi).standardize("submitted name", taxon_id_hint=12345)

    assert ncbi.name_requests == []
    assert record.taxon_id == 12345
    assert record.scientific_name == "Iris japonica"


def test_ete_resolver_keeps_current_infraspecific_taxon_name():
    class SubspeciesNCBI(FakeNCBITaxa):
        def get_lineage(self, taxid: int) -> list[int]:
            self.lineage_requests.append(taxid)
            return [33090, 35493, 3398, 73496, 58920, 26379, 12345, 54321]

        def get_rank(self, taxids: list[int]) -> dict[int, str]:
            self.rank_requests.append(taxids)
            ranks = super().get_rank(taxids)
            ranks[54321] = "subspecies"
            return ranks

        def get_taxid_translator(self, taxids: list[int]) -> dict[int, str]:
            self.taxid_requests.append(taxids)
            names = super().get_taxid_translator(taxids)
            names[54321] = "Iris japonica subsp. formosana"
            return names

    record = ETETaxonomyResolver(ncbi=SubspeciesNCBI()).standardize(
        "submitted name",
        taxon_id_hint=54321,
    )

    assert record.taxon_id == 54321
    assert record.scientific_name == "Iris japonica subsp. formosana"
    assert record.species == "japonica"
    assert record.infraspecific_rank == "subspecies"


def test_ete_resolver_raises_when_name_cannot_be_resolved():
    class UnresolvedNCBI(FakeNCBITaxa):
        def get_name_translator(self, names: list[str]) -> dict[str, list[int]]:
            return {}

    with pytest.raises(TaxonomyError, match="ETE could not resolve taxon name"):
        ETETaxonomyResolver(ncbi=UnresolvedNCBI()).standardize("Missing species")
