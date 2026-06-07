from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from barcode_kit import cli
from barcode_kit.config import AppConfig, TreeShrinkConfig
from barcode_kit.models import BuildReportEntry, Marker, TaxonExclusion, TaxonQuery


def test_build_creates_missing_config_file_before_running(tmp_path: Path, monkeypatch):
    def fake_build_dataset(
        config: AppConfig,
        storage: Any,
        query: TaxonQuery,
        marker: Marker,
        outdir: Path,
        **kwargs: Any,
    ) -> list[BuildReportEntry]:
        assert config.collectors.batch_size == 500
        assert query == TaxonQuery("genus", "Iris")
        assert marker is Marker.ITS
        assert kwargs == {}
        return []

    config_path = tmp_path / "config.toml"
    data_dir = tmp_path / "data"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    monkeypatch.setattr(cli.config_module, "DEFAULT_DATA_DIR", data_dir)
    monkeypatch.setattr(cli, "Storage", lambda path: object())
    monkeypatch.setattr(cli, "build_dataset", fake_build_dataset)

    result = CliRunner().invoke(
        cli.app,
        ["build", "--genus", "Iris", "--marker", "its", "--outdir", str(tmp_path / "out")],
    )

    assert result.exit_code == 0, result.output
    assert config_path.exists()
    assert "[build.blast_rescue]" in config_path.read_text(encoding="utf-8")


def test_build_rejects_removed_its_extraction_mode_option(tmp_path: Path, monkeypatch):
    config_path = tmp_path / "config.toml"
    data_dir = tmp_path / "data"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    monkeypatch.setattr(cli.config_module, "DEFAULT_DATA_DIR", data_dir)

    result = CliRunner().invoke(
        cli.app,
        [
            "build",
            "--genus",
            "Iris",
            "--marker",
            "its",
            "--its-extraction-mode",
            "annotation",
        ],
    )

    assert result.exit_code != 0
    assert "No such option" in result.output


def test_build_accepts_lineage_constraints_for_taxon_query(tmp_path: Path, monkeypatch):
    def fake_build_dataset(
        config: AppConfig,
        storage: Any,
        query: TaxonQuery,
        marker: Marker,
        outdir: Path,
        **kwargs: Any,
    ) -> list[BuildReportEntry]:
        assert query.rank == "genus"
        assert query.name == "Iris"
        assert [(item.rank, item.name) for item in query.constraints] == [
            ("kingdom", "Viridiplantae"),
        ]
        assert marker is Marker.RBCL
        return []

    config_path = tmp_path / "config.toml"
    data_dir = tmp_path / "data"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    monkeypatch.setattr(cli.config_module, "DEFAULT_DATA_DIR", data_dir)
    monkeypatch.setattr(cli, "Storage", lambda path: object())
    monkeypatch.setattr(cli, "build_dataset", fake_build_dataset)

    result = CliRunner().invoke(
        cli.app,
        [
            "build",
            "--kindom",
            "Viridiplantae",
            "--genus",
            "Iris",
            "--marker",
            "rbcl",
            "--outdir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 0, result.output


@pytest.mark.parametrize(
    "option_args",
    [
        ["--exclude", "hybrid"],
        ["--min-length", "500"],
        ["--max-ambiguous", "0.05"],
        ["--tree-shrink-qc"],
    ],
)
def test_build_rejects_moved_qc_options(option_args):
    result = CliRunner().invoke(
        cli.app,
        [
            "build",
            "--genus",
            "Iris",
            "--marker",
            "rbcl",
            *option_args,
        ],
    )

    assert result.exit_code != 0
    assert "No such option" in result.output


def test_build_rejects_removed_legacy_exclude_options():
    for option in ("--exclude-hybrid", "--exclude-uncertain"):
        result = CliRunner().invoke(
            cli.app,
            ["build", "--genus", "Iris", "--marker", "rbcl", option],
        )

        assert result.exit_code != 0
        assert "No such option" in result.output


@pytest.mark.parametrize(
    ("option_args", "expected"),
    [
        (
            ["--exclude", "hybrid", "--exclude", "variety"],
            {
                "exclude": {TaxonExclusion.HYBRID, TaxonExclusion.VARIETY},
                "min_length": None,
                "max_ambiguous_content": None,
                "enable_tree_shrink_qc": False,
            },
        ),
        (
            ["--min-length", "500"],
            {
                "exclude": set(),
                "min_length": 500,
                "max_ambiguous_content": None,
                "enable_tree_shrink_qc": False,
            },
        ),
        (
            ["--max-ambiguous", "0.05"],
            {
                "exclude": set(),
                "min_length": None,
                "max_ambiguous_content": 0.05,
                "enable_tree_shrink_qc": False,
            },
        ),
        (
            ["--tree-shrink-qc"],
            {
                "exclude": set(),
                "min_length": None,
                "max_ambiguous_content": None,
                "enable_tree_shrink_qc": True,
            },
        ),
    ],
)
def test_qc_runs_each_check_independently(
    tmp_path: Path,
    monkeypatch,
    option_args,
    expected,
):
    qc_calls = []

    def fake_run_qc(
        dataset: Path,
        *,
        exclude: set[TaxonExclusion],
        min_length: int | None,
        max_ambiguous_content: float | None,
        enable_tree_shrink_qc: bool,
        tree_shrink_config: TreeShrinkConfig,
    ) -> dict[str, Any]:
        qc_calls.append(
            {
                "dataset": dataset,
                "exclude": exclude,
                "min_length": min_length,
                "max_ambiguous_content": max_ambiguous_content,
                "enable_tree_shrink_qc": enable_tree_shrink_qc,
                "tree_shrink_config": tree_shrink_config,
            }
        )
        return {"records": []}

    config_path = tmp_path / "config.toml"
    data_dir = tmp_path / "data"
    config_path.write_text(
        f"""
[paths]
data_dir = "{data_dir}"

[build.tree_shrink_qc]
quantile = 0.01
bootstrap = 1000
max_removed = 6
""",
        encoding="utf-8",
    )
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    monkeypatch.setattr(cli.validation_module, "run_qc", fake_run_qc, raising=False)

    result = CliRunner().invoke(
        cli.app,
        [
            "qc",
            "--dataset",
            str(tmp_path / "dataset"),
            *option_args,
        ],
    )

    assert result.exit_code == 0, result.output
    assert qc_calls == [
        {
            "dataset": tmp_path / "dataset",
            **expected,
            "tree_shrink_config": TreeShrinkConfig(
                quantile=0.01,
                bootstrap=1000,
                max_removed=6,
            ),
        }
    ]


def test_qc_requires_at_least_one_check(tmp_path: Path, monkeypatch):
    config_path = tmp_path / "config.toml"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))

    result = CliRunner().invoke(
        cli.app,
        [
            "qc",
            "--dataset",
            str(tmp_path / "dataset"),
        ],
    )

    assert result.exit_code != 0
    assert "select at least one QC option" in result.output


def test_qc_rejects_treeshrink_qc_threads_cli_option(tmp_path: Path, monkeypatch):
    config_path = tmp_path / "config.toml"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))

    result = CliRunner().invoke(
        cli.app,
        [
            "qc",
            "--dataset",
            str(tmp_path / "dataset"),
            "--tree-shrink-qc-threads",
            "8",
        ],
    )

    assert result.exit_code != 0
    assert "No such option: --tree-shrink-qc-threads" in result.output
