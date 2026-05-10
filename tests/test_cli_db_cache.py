from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from barcode_kit.cli import app
from barcode_kit.config import AppConfig, write_config
from barcode_kit.models import GenBankCacheRecord, TaxonomyRecord
from barcode_kit.storage import Storage


def test_db_remove_by_accession_deletes_database_row_and_file(tmp_path, monkeypatch):
    config = _write_test_config(tmp_path, monkeypatch)
    storage = Storage(config.database_path)
    _seed_record(storage, "PP476489", 4, 12345, "Iris japonica")
    _seed_record(storage, "AB000001", 2, 67890, "Rosa chinensis")
    _cache_file(config, "PP476489.4").write_text("cached iris", encoding="utf-8")
    _cache_file(config, "AB000001.2").write_text("cached rosa", encoding="utf-8")

    result = CliRunner().invoke(app, ["db", "remove", "--accession", "PP476489.4", "--yes"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["database_records_removed"] == 1
    assert payload["files_removed"] == 1
    assert payload["taxonomy_removed"] == 1
    assert not _cache_file(config, "PP476489.4").exists()
    assert _cache_file(config, "AB000001.2").exists()
    assert storage.counts() == {"taxonomy": 1, "genbank_cache": 1}


def test_db_info_rank_lists_species_for_genus(tmp_path, monkeypatch):
    config = _write_test_config(tmp_path, monkeypatch)
    storage = Storage(config.database_path)
    _seed_record(
        storage,
        "PX743804",
        1,
        12345,
        "Aspidistra elatior",
        family="Asparagaceae",
        genus="Aspidistra",
    )
    _seed_record(
        storage,
        "PX743803",
        1,
        67890,
        "Aspidistra typica",
        family="Asparagaceae",
        genus="Aspidistra",
    )
    _seed_record(
        storage,
        "AB000001",
        2,
        24680,
        "Rosa chinensis",
        family="Rosaceae",
        genus="Rosa",
    )

    result = CliRunner().invoke(app, ["db", "info", "--rank", "species", "--genus", "Aspidistra"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["rank"] == "species"
    assert payload["query"] == {"rank": "genus", "name": "Aspidistra"}
    assert payload["taxa"] == [
        {
            "name": "Aspidistra elatior",
            "records": 1,
            "markers": {
                "its": 0,
                "its2": 0,
                "matk": 0,
                "rbcl": 0,
            },
        },
        {
            "name": "Aspidistra typica",
            "records": 1,
            "markers": {
                "its": 0,
                "its2": 0,
                "matk": 0,
                "rbcl": 0,
            },
        },
    ]


def test_db_info_rank_genus_requires_family_filter_when_filtering(tmp_path, monkeypatch):
    _write_test_config(tmp_path, monkeypatch)

    result = CliRunner().invoke(app, ["db", "info", "--rank", "genus", "--genus", "Aspidistra"])

    assert result.exit_code == 1
    assert "--rank genus can only be filtered by --family" in result.output


def test_db_remove_by_genus_leaves_orphan_genbank_file_for_prune(
    tmp_path,
    monkeypatch,
    genbank_text,
):
    config = _write_test_config(tmp_path, monkeypatch)
    storage = Storage(config.database_path)
    storage.initialize()
    _cache_file(config, "PX743804.1").write_text(
        genbank_text(accession="PX743804", version=1, organism="Aspidistra elatior"),
        encoding="utf-8",
    )
    _cache_file(config, "AB000001.2").write_text(
        genbank_text(accession="AB000001", version=2, organism="Rosa chinensis"),
        encoding="utf-8",
    )

    result = CliRunner().invoke(app, ["db", "remove", "--genus", "Aspidistra", "--yes"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["database_records_removed"] == 0
    assert payload["files_removed"] == 0
    assert payload["taxonomy_removed"] == 0
    assert _cache_file(config, "PX743804.1").exists()
    assert _cache_file(config, "AB000001.2").exists()


def test_db_clear_deletes_all_cache_rows_taxonomy_and_files(tmp_path, monkeypatch):
    config = _write_test_config(tmp_path, monkeypatch)
    storage = Storage(config.database_path)
    _seed_record(storage, "PP476489", 4, 12345, "Iris japonica")
    _seed_record(storage, "AB000001", 2, 67890, "Rosa chinensis")
    _cache_file(config, "PP476489.4").write_text("cached iris", encoding="utf-8")
    _cache_file(config, "AB000001.2").write_text("cached rosa", encoding="utf-8")

    result = CliRunner().invoke(app, ["db", "clear", "--yes"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["database_records_removed"] == 2
    assert payload["taxonomy_removed"] == 2
    assert payload["files_removed"] == 2
    assert storage.counts() == {"taxonomy": 0, "genbank_cache": 0}
    assert list(config.genbank_cache_dir.glob("*.gb")) == []


def test_db_prune_removes_missing_rows_unreferenced_files_and_orphan_taxonomy(
    tmp_path,
    monkeypatch,
):
    config = _write_test_config(tmp_path, monkeypatch)
    storage = Storage(config.database_path)
    _seed_record(storage, "PP476489", 4, 12345, "Iris japonica")
    _seed_record(storage, "AB000001", 2, 67890, "Rosa chinensis")
    _cache_file(config, "AB000001.2").write_text("referenced", encoding="utf-8")
    _cache_file(config, "OLD000001.1").write_text("unreferenced", encoding="utf-8")

    result = CliRunner().invoke(app, ["db", "prune", "--yes"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["database_records_removed"] == 1
    assert payload["files_removed"] == 1
    assert payload["taxonomy_removed"] == 1
    assert storage.counts() == {"taxonomy": 1, "genbank_cache": 1}
    assert not _cache_file(config, "PP476489.4").exists()
    assert _cache_file(config, "AB000001.2").exists()
    assert not _cache_file(config, "OLD000001.1").exists()


def test_db_prune_removes_orphan_files_when_database_is_empty(
    tmp_path,
    monkeypatch,
    genbank_text,
):
    config = _write_test_config(tmp_path, monkeypatch)
    Storage(config.database_path).initialize()
    _cache_file(config, "PX743804.1").write_text(
        genbank_text(accession="PX743804", version=1, organism="Aspidistra elatior"),
        encoding="utf-8",
    )
    _cache_file(config, "PX743803.1").write_text(
        genbank_text(accession="PX743803", version=1, organism="Aspidistra typica"),
        encoding="utf-8",
    )

    result = CliRunner().invoke(app, ["db", "prune", "--yes"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["database_records_removed"] == 0
    assert payload["files_removed"] == 2
    assert payload["taxonomy_removed"] == 0
    assert list(config.genbank_cache_dir.glob("*.gb")) == []


def _write_test_config(tmp_path: Path, monkeypatch) -> AppConfig:
    config_path = tmp_path / "config.toml"
    config = AppConfig(
        data_dir=tmp_path / "data",
        batch_size=500,
        download_workers=1,
        timeout=30,
        retry_attempts=3,
        genbank_email="test@example.com",
    )
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    write_config(config, config_path)
    return config


def _seed_record(
    storage: Storage,
    accession_root: str,
    version: int,
    taxon_id: int,
    scientific_name: str,
    family: str | None = None,
    genus: str | None = None,
) -> None:
    storage.initialize()
    storage.upsert_taxonomy(
        TaxonomyRecord(
            taxon_id=taxon_id,
            scientific_name=scientific_name,
            family=family,
            genus=genus,
        )
    )
    storage.upsert_genbank_cache(
        GenBankCacheRecord(
            accession_root=accession_root,
            version=version,
            accession_version=f"{accession_root}.{version}",
            taxon_id=taxon_id,
        )
    )


def _cache_file(config: AppConfig, accession_version: str) -> Path:
    config.genbank_cache_dir.mkdir(parents=True, exist_ok=True)
    return config.genbank_cache_dir / f"{accession_version}.gb"
