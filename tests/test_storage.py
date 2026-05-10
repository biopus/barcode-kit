from __future__ import annotations

from barcode_kit.models import GenBankCacheRecord, TaxonQuery, TaxonomyRecord
from barcode_kit.storage import Storage


def test_get_cached_versions_by_roots_returns_only_cached_versions(tmp_path):
    storage = Storage(tmp_path / "cache.db")
    storage.initialize()
    storage.upsert_taxonomy(TaxonomyRecord(taxon_id=12345, scientific_name="Iris japonica"))
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="PP476489",
            version=4,
            accession_version="PP476489.4",
            taxon_id=12345,
            has_rbcl=True,
        )
    )
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="AB000001",
            version=2,
            accession_version="AB000001.2",
            taxon_id=12345,
            has_matk=True,
        )
    )

    versions = storage.get_cached_versions_by_roots(
        ["PP476489", "missing", "PP476489", "AB000001"]
    )

    assert versions == {"PP476489": 4, "AB000001": 2}


def test_cache_records_filters_by_taxon_and_accession(tmp_path):
    storage = Storage(tmp_path / "cache.db")
    storage.initialize()
    storage.upsert_taxonomy(
        TaxonomyRecord(
            taxon_id=12345,
            scientific_name="Iris japonica",
            family="Iridaceae",
            genus="Iris",
            species="Iris japonica",
        )
    )
    storage.upsert_taxonomy(
        TaxonomyRecord(
            taxon_id=67890,
            scientific_name="Rosa chinensis",
            family="Rosaceae",
            genus="Rosa",
            species="Rosa chinensis",
        )
    )
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="PP476489",
            version=4,
            accession_version="PP476489.4",
            taxon_id=12345,
            has_rbcl=True,
        )
    )
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="AB000001",
            version=2,
            accession_version="AB000001.2",
            taxon_id=67890,
            has_matk=True,
        )
    )

    iris_records = storage.cache_records(TaxonQuery("genus", "Iris"))
    accession_records = storage.cache_records(accession="AB000001.2")

    assert [record.accession_version for record in iris_records] == ["PP476489.4"]
    assert [record.accession_root for record in accession_records] == ["AB000001"]


def test_taxon_summaries_lists_genus_with_counts_and_marker_coverage(tmp_path):
    storage = Storage(tmp_path / "cache.db")
    storage.initialize()
    storage.upsert_taxonomy(
        TaxonomyRecord(
            taxon_id=12345,
            scientific_name="Iris japonica",
            family="Iridaceae",
            genus="Iris",
        )
    )
    storage.upsert_taxonomy(
        TaxonomyRecord(
            taxon_id=67890,
            scientific_name="Iris tectorum",
            family="Iridaceae",
            genus="Iris",
        )
    )
    storage.upsert_taxonomy(
        TaxonomyRecord(
            taxon_id=24680,
            scientific_name="Aspidistra elatior",
            family="Asparagaceae",
            genus="Aspidistra",
        )
    )
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="PP476489",
            version=4,
            accession_version="PP476489.4",
            taxon_id=12345,
            has_rbcl=True,
        )
    )
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="PP476490",
            version=1,
            accession_version="PP476490.1",
            taxon_id=67890,
            has_matk=True,
        )
    )
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="PX743804",
            version=1,
            accession_version="PX743804.1",
            taxon_id=24680,
            has_its=True,
        )
    )

    summaries = storage.taxon_summaries("genus", TaxonQuery("family", "Iridaceae"))

    assert summaries == [
        {
            "name": "Iris",
            "records": 2,
            "markers": {
                "its": 0,
                "its2": 0,
                "matk": 1,
                "rbcl": 1,
            },
        }
    ]


def test_delete_cache_records_removes_rows_and_orphan_taxonomy(tmp_path):
    storage = Storage(tmp_path / "cache.db")
    storage.initialize()
    storage.upsert_taxonomy(TaxonomyRecord(taxon_id=12345, scientific_name="Iris japonica"))
    storage.upsert_taxonomy(TaxonomyRecord(taxon_id=67890, scientific_name="Rosa chinensis"))
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="PP476489",
            version=4,
            accession_version="PP476489.4",
            taxon_id=12345,
        )
    )
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="AB000001",
            version=2,
            accession_version="AB000001.2",
            taxon_id=67890,
        )
    )

    deleted_cache = storage.delete_cache_records(["PP476489"])
    deleted_taxonomy = storage.delete_orphan_taxonomy()

    assert deleted_cache == 1
    assert deleted_taxonomy == 1
    assert storage.counts() == {"taxonomy": 1, "genbank_cache": 1}
    assert storage.cache_records(accession="PP476489.4") == []


def test_clear_cache_removes_all_cache_and_taxonomy(tmp_path):
    storage = Storage(tmp_path / "cache.db")
    storage.initialize()
    storage.upsert_taxonomy(TaxonomyRecord(taxon_id=12345, scientific_name="Iris japonica"))
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root="PP476489",
            version=4,
            accession_version="PP476489.4",
            taxon_id=12345,
        )
    )

    result = storage.clear_cache()

    assert result == {"genbank_cache": 1, "taxonomy": 1}
    assert storage.counts() == {"taxonomy": 0, "genbank_cache": 0}
