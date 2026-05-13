from __future__ import annotations

from pathlib import Path
from typing import Any

from typer.testing import CliRunner

from barcode_kit import cli
from barcode_kit.config import AppConfig
from barcode_kit.models import BuildReportEntry, Marker, TaxonQuery
from barcode_kit.phylogeny import TreeShrinkQcConfig


def test_build_creates_missing_config_file_before_running(tmp_path: Path, monkeypatch):
    def fake_build_dataset(
        config: AppConfig,
        storage: Any,
        query: TaxonQuery,
        marker: Marker,
        outdir: Path,
        **kwargs: Any,
    ) -> list[BuildReportEntry]:
        assert config.batch_size == 500
        assert query == TaxonQuery("genus", "Iris")
        assert marker is Marker.ITS
        return []

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
    monkeypatch.setattr(cli, "build_dataset", fake_build_dataset)

    result = CliRunner().invoke(
        cli.app,
        ["build", "--genus", "Iris", "--marker", "its", "--outdir", str(tmp_path / "out")],
    )

    assert result.exit_code == 0, result.output
    assert config_path.exists()
    assert "[build.blast_rescue]" in config_path.read_text(encoding="utf-8")


def test_build_passes_treeshrink_qc_options(tmp_path: Path, monkeypatch):
    def fake_build_dataset(
        config: AppConfig,
        storage: Any,
        query: TaxonQuery,
        marker: Marker,
        outdir: Path,
        **kwargs: Any,
    ) -> list[BuildReportEntry]:
        assert kwargs["tree_shrink_qc"] == TreeShrinkQcConfig(threads=8, quantile=0.01)
        return []

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
    monkeypatch.setattr(cli, "build_dataset", fake_build_dataset)

    result = CliRunner().invoke(
        cli.app,
        [
            "build",
            "--genus",
            "Iris",
            "--marker",
            "rbcl",
            "--outdir",
            str(tmp_path / "out"),
            "--tree-shrink-qc",
            "--tree-shrink-qc-threads",
            "8",
            "--tree-shrink-qc-quantile",
            "0.01",
        ],
    )

    assert result.exit_code == 0, result.output
