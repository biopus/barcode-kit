from __future__ import annotations

from typer.testing import CliRunner

from barcode_kit import config as config_module
from barcode_kit.cli import app


def test_cli_does_not_expose_config_command():
    runner = CliRunner()

    help_result = runner.invoke(app, ["--help"])
    config_result = runner.invoke(app, ["config", "list"])

    assert help_result.exit_code == 0, help_result.output
    assert "config" not in help_result.output
    assert config_result.exit_code != 0
    assert "No such command" in config_result.output


def test_load_or_create_config_writes_blast_rescue_defaults(tmp_path, monkeypatch):
    config_path = tmp_path / "config.toml"
    data_dir = tmp_path / "data"
    monkeypatch.setattr(config_module, "DEFAULT_DATA_DIR", data_dir)

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


def test_load_config_reads_partial_file_with_dataclass_defaults(tmp_path):
    config_path = tmp_path / "config.toml"
    config_path.write_text(
        "\n".join(
            [
                "[collectors.genbank]",
                'email = "user@example.com"',
                "",
                "[build.blast_rescue.its]",
                "min_identity = 0.92",
                "",
            ]
        ),
        encoding="utf-8",
    )

    config = config_module.load_config(config_path)

    assert config.collectors.genbank_email == "user@example.com"
    assert config.collectors.batch_size == 500
    assert config.collectors.download_workers == 8
    assert config.blast_rescue.its.min_identity == 0.92
    assert config.blast_rescue.its.min_subject_coverage == 0.85
    assert config.blast_rescue.its.min_query_length_ratio == 0.85
    assert config.blast_rescue.its.max_query_length_ratio == 1.20
    assert not hasattr(config_module, "_deep_merge")
    assert not hasattr(config_module, "DEFAULT_CONFIG")
