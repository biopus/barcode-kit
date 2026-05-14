from __future__ import annotations

from typer.testing import CliRunner

from barcode_kit import config as config_module
from barcode_kit.cli import app


def test_config_set_and_list(tmp_path, monkeypatch):
    config_path = tmp_path / "config.toml"
    data_dir = tmp_path / "data"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    runner = CliRunner()

    result = runner.invoke(app, ["config", "set", "paths.data_dir", str(data_dir)])
    assert result.exit_code == 0, result.output

    result = runner.invoke(app, ["config", "set", "collectors.genbank.email", "user@example.com"])
    assert result.exit_code == 0, result.output

    result = runner.invoke(app, ["config", "list"])
    assert result.exit_code == 0, result.output
    assert "user@example.com" in result.output
    assert str(data_dir) in result.output


def test_load_or_create_config_writes_blast_rescue_defaults(tmp_path, monkeypatch):
    config_path = tmp_path / "config.toml"
    data_dir = tmp_path / "data"
    monkeypatch.setattr(
        config_module,
        "DEFAULT_CONFIG",
        {
            **config_module.DEFAULT_CONFIG,
            "paths": {"data_dir": str(data_dir)},
        },
    )

    config = config_module.load_or_create_config(config_path)

    assert config_path.exists()
    assert config.data_dir == data_dir
    assert config.blast_rescue.its.min_identity == 0.80
    assert config.blast_rescue.its2.min_subject_coverage == 0.90
    assert config.blast_rescue.word_size == 11
    assert config.blast_rescue.evalue == 1e-3
    assert config.itsxrust.inc_e == 0.01
    assert config.itsxrust.min_anchor_score == 8
    assert config.itsxrust.max_per_anchor == 20
    assert config.itsxrust.max_anchor_evalue == 0.01
    assert config.tree_shrink_qc.quantile == 0.1
    assert config.tree_shrink_qc.bootstrap == 0
    assert config.tree_shrink_qc.max_removed is None
    contents = config_path.read_text(encoding="utf-8")
    assert "[build.blast_rescue]" in contents
    assert "[build.blast_rescue.its]" in contents
    assert "[build.blast_rescue.its2]" in contents
    assert "[build.itsxrust]" in contents
    assert "[build.tree_shrink_qc]" in contents
    assert 'max_removed = "auto-select"' in contents
    assert "makeblastdb_command" not in contents
    assert "blastn_command" not in contents
    assert "blastn_outfmt" not in contents


def test_config_set_updates_blast_rescue_values(tmp_path, monkeypatch):
    config_path = tmp_path / "config.toml"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["config", "set", "build.blast_rescue.its.min_identity", "0.92"],
    )

    assert result.exit_code == 0, result.output
    config = config_module.load_config(config_path)
    assert config.blast_rescue.its.min_identity == 0.92


def test_config_set_updates_blast_and_itsxrust_values(tmp_path, monkeypatch):
    config_path = tmp_path / "config.toml"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    runner = CliRunner()

    blast_result = runner.invoke(
        app,
        ["config", "set", "build.blast_rescue.word_size", "15"],
    )
    itsxrust_result = runner.invoke(
        app,
        ["config", "set", "build.itsxrust.max_per_anchor", "30"],
    )

    assert blast_result.exit_code == 0, blast_result.output
    assert itsxrust_result.exit_code == 0, itsxrust_result.output
    config = config_module.load_config(config_path)
    assert config.blast_rescue.word_size == 15
    assert config.itsxrust.max_per_anchor == 30


def test_config_set_updates_tree_shrink_qc_values(tmp_path, monkeypatch):
    config_path = tmp_path / "config.toml"
    monkeypatch.setenv("BARCODE_KIT_CONFIG", str(config_path))
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["config", "set", "build.tree_shrink_qc.quantile", "0.01"],
    )
    bootstrap_result = runner.invoke(
        app,
        ["config", "set", "build.tree_shrink_qc.bootstrap", "1000"],
    )
    max_removed_result = runner.invoke(
        app,
        ["config", "set", "build.tree_shrink_qc.max_removed", "6"],
    )
    auto_select_result = runner.invoke(
        app,
        ["config", "set", "build.tree_shrink_qc.max_removed", "auto-select"],
    )

    assert result.exit_code == 0, result.output
    assert bootstrap_result.exit_code == 0, bootstrap_result.output
    assert max_removed_result.exit_code == 0, max_removed_result.output
    assert auto_select_result.exit_code == 0, auto_select_result.output
    config = config_module.load_config(config_path)
    assert config.tree_shrink_qc.quantile == 0.01
    assert config.tree_shrink_qc.bootstrap == 1000
    assert config.tree_shrink_qc.max_removed is None
    assert 'max_removed = "auto-select"' in config_path.read_text(encoding="utf-8")
