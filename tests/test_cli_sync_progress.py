from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from typer.testing import CliRunner

from barcode_kit import cli
from barcode_kit.config import AppConfig
from barcode_kit.genbank import SyncResult
from barcode_kit.models import Marker, TaxonQuery


def test_sync_writes_live_download_rate_to_stderr_and_final_json_to_stdout(
    tmp_path: Path,
    monkeypatch,
):
    class FakeSyncService:
        def __init__(self, config: AppConfig, storage: Any, taxonomy_resolver: Any):
            pass

        def sync(
            self,
            query: TaxonQuery,
            marker: Marker,
            progress: Any | None = None,
        ) -> SyncResult:
            assert query == TaxonQuery("genus", "Iris")
            assert marker is Marker.RBCL
            assert progress is not None
            progress.start_download(total_records=1)
            progress.record_downloaded_bytes(2048)
            progress.record_downloaded_record("PP476489.4")
            return SyncResult(
                query="Iris[Organism] AND rbcl",
                remote_count=1,
                downloaded=1,
                reused_local=0,
                ingested=1,
                skipped=0,
                updated=0,
                failed=[],
            )

        def close(self) -> None:
            pass

    config = AppConfig(
        data_dir=tmp_path / "data",
        batch_size=500,
        download_workers=1,
        timeout=30,
        retry_attempts=3,
        genbank_email="test@example.com",
    )
    monkeypatch.setattr(cli.config_module, "load_or_create_config", lambda: config)
    monkeypatch.setattr(cli.config_module, "ensure_app_dirs", lambda config: None)
    monkeypatch.setattr(cli, "Storage", lambda path: object())
    monkeypatch.setattr(cli, "ETETaxonomyResolver", lambda: object())
    monkeypatch.setattr(cli, "SyncService", FakeSyncService)

    result = CliRunner().invoke(
        cli.app,
        ["sync", "--genus", "Iris", "--marker", "rbcl"],
    )

    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout) == asdict(
        SyncResult(
            query="Iris[Organism] AND rbcl",
            remote_count=1,
            downloaded=1,
            reused_local=0,
            ingested=1,
            skipped=0,
            updated=0,
            failed=[],
        )
    )
    assert "downloaded=1/1" in result.stderr
    assert "speed=" in result.stderr
    assert "avg=" in result.stderr


def test_sync_creates_missing_config_file_before_running(tmp_path: Path, monkeypatch):
    class FakeSyncService:
        def __init__(self, config: AppConfig, storage: Any, taxonomy_resolver: Any):
            assert config.batch_size == 500

        def sync(
            self,
            query: TaxonQuery,
            marker: Marker,
            progress: Any | None = None,
        ) -> SyncResult:
            return SyncResult(
                query="Iris[Organism] AND rbcl",
                remote_count=0,
                downloaded=0,
                reused_local=0,
                ingested=0,
                skipped=0,
                updated=0,
                failed=[],
            )

        def close(self) -> None:
            pass

    config_path = tmp_path / "config.toml"
    data_dir = tmp_path / "data"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    monkeypatch.setattr(
        cli.config_module,
        "DEFAULT_CONFIG",
        {
            **cli.config_module.DEFAULT_CONFIG,
            "paths": {"data_dir": str(data_dir)},
        },
    )
    monkeypatch.setattr(cli, "Storage", lambda path: object())
    monkeypatch.setattr(cli, "ETETaxonomyResolver", lambda: object())
    monkeypatch.setattr(cli, "SyncService", FakeSyncService)

    result = CliRunner().invoke(
        cli.app,
        ["sync", "--genus", "Iris", "--marker", "rbcl"],
    )

    assert result.exit_code == 0, result.output
    assert config_path.exists()
    assert "[build.blast_rescue]" in config_path.read_text(encoding="utf-8")


def test_sync_rejects_trnl_trnf_marker(tmp_path: Path, monkeypatch):
    class FailingSyncService:
        def __init__(self, config: AppConfig, storage: Any, taxonomy_resolver: Any):
            raise AssertionError("trnl-trnf should be rejected before sync starts")

    config = AppConfig(
        data_dir=tmp_path / "data",
        batch_size=500,
        download_workers=1,
        timeout=30,
        retry_attempts=3,
        genbank_email="test@example.com",
    )
    monkeypatch.setattr(cli.config_module, "load_config", lambda: config)
    monkeypatch.setattr(cli.config_module, "ensure_app_dirs", lambda config: None)
    monkeypatch.setattr(cli, "Storage", lambda path: object())
    monkeypatch.setattr(cli, "ETETaxonomyResolver", lambda: object())
    monkeypatch.setattr(cli, "SyncService", FailingSyncService)

    result = CliRunner().invoke(
        cli.app,
        ["sync", "--genus", "Iris", "--marker", "trnl-trnf"],
    )

    assert result.exit_code == 2
    assert "invalid" in result.output.lower()
