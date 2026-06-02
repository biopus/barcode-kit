from __future__ import annotations

import os
import tomllib
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any

from barcode_kit.exceptions import ConfigError


__all__ = [
    "AppConfig",
    "BlastRescueConfig",
    "BlastRescueMarkerConfig",
    "CollectorConfig",
    "DEFAULT_DATA_DIR",
    "ItsxrustConfig",
    "TreeShrinkConfig",
    "config_as_dict",
    "ensure_app_dirs",
    "load_config",
    "load_or_create_config",
    "write_config",
]


DEFAULT_DATA_DIR = Path.home() / ".barcode-kit"


@dataclass(frozen=True)
class BlastRescueMarkerConfig:
    min_subject_coverage: float
    min_identity: float
    min_query_length_ratio: float
    max_query_length_ratio: float


@dataclass(frozen=True)
class BlastRescueConfig:
    blastn_dust: str = "no"
    word_size: int = 11
    evalue: float = 1e-3
    endpoint_margin_bases: int = 15
    endpoint_margin_fraction: float = 0.05
    ambiguous_bitscore_ratio: float = 0.95
    ambiguous_overlap_ratio: float = 0.80
    its: BlastRescueMarkerConfig = field(
        default_factory=lambda: BlastRescueMarkerConfig(
            min_subject_coverage=0.85,
            min_identity=0.80,
            min_query_length_ratio=0.85,
            max_query_length_ratio=1.20,
        )
    )
    its2: BlastRescueMarkerConfig = field(
        default_factory=lambda: BlastRescueMarkerConfig(
            min_subject_coverage=0.90,
            min_identity=0.85,
            min_query_length_ratio=0.90,
            max_query_length_ratio=1.15,
        )
    )


@dataclass(frozen=True)
class ItsxrustConfig:
    inc_e: float = 0.01
    min_anchor_score: int = 8
    max_per_anchor: int = 20
    max_anchor_evalue: float = 0.01


@dataclass(frozen=True)
class TreeShrinkConfig:
    quantile: float = 0.1
    bootstrap: int = 0
    max_removed: int | None = None


@dataclass(frozen=True)
class CollectorConfig:
    batch_size: int = 500
    download_workers: int = 8
    timeout: float = 30
    retry_attempts: int = 3
    genbank_email: str = ""
    genbank_api_key: str | None = None


@dataclass(frozen=True)
class AppConfig:
    data_dir: Path = field(default_factory=lambda: DEFAULT_DATA_DIR)
    collectors: CollectorConfig = field(default_factory=CollectorConfig)
    blast_rescue: BlastRescueConfig = field(default_factory=BlastRescueConfig)
    itsxrust: ItsxrustConfig = field(default_factory=ItsxrustConfig)
    tree_shrink_qc: TreeShrinkConfig = field(default_factory=TreeShrinkConfig)

    @property
    def database_path(self) -> Path:
        return self.data_dir / "database.db"

    @property
    def config_path(self) -> Path:
        return self.data_dir / "config.toml"

    @property
    def genbank_cache_dir(self) -> Path:
        return self.data_dir / "cache" / "genbank"

    @property
    def logs_dir(self) -> Path:
        return self.data_dir / "logs"


def load_config(path: Path | None = None) -> AppConfig:
    path = path or Path(os.environ.get("BARCODE_KIT_CONFIG", DEFAULT_DATA_DIR / "config.toml"))
    if not path.exists():
        return AppConfig()
    return _parse_config(_read_toml(path))


def load_or_create_config(path: Path | None = None) -> AppConfig:
    path = path or Path(os.environ.get("BARCODE_KIT_CONFIG", DEFAULT_DATA_DIR / "config.toml"))
    config = load_config(path)
    ensure_app_dirs(config)
    if not path.exists():
        write_config(config, path)
    return config


def ensure_app_dirs(config: AppConfig) -> None:
    config.data_dir.mkdir(parents=True, exist_ok=True)
    config.genbank_cache_dir.mkdir(parents=True, exist_ok=True)
    config.logs_dir.mkdir(parents=True, exist_ok=True)


def write_config(config: AppConfig, path: Path | None = None) -> None:
    path = path or config.config_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_toml_dump(config_as_dict(config)), encoding="utf-8")


def config_as_dict(config: AppConfig) -> dict[str, Any]:
    return {
        "paths": {"data_dir": str(config.data_dir)},
        "collectors": {
            "batch_size": config.collectors.batch_size,
            "download_workers": config.collectors.download_workers,
            "timeout": config.collectors.timeout,
            "retry_attempts": config.collectors.retry_attempts,
            "genbank": {
                "email": config.collectors.genbank_email,
                "api_key": config.collectors.genbank_api_key or "",
            },
        },
        "build": {
            "itsxrust": {
                "inc_e": config.itsxrust.inc_e,
                "min_anchor_score": config.itsxrust.min_anchor_score,
                "max_per_anchor": config.itsxrust.max_per_anchor,
                "max_anchor_evalue": config.itsxrust.max_anchor_evalue,
            },
            "blast_rescue": {
                "blastn_dust": config.blast_rescue.blastn_dust,
                "word_size": config.blast_rescue.word_size,
                "evalue": config.blast_rescue.evalue,
                "endpoint_margin_bases": config.blast_rescue.endpoint_margin_bases,
                "endpoint_margin_fraction": config.blast_rescue.endpoint_margin_fraction,
                "ambiguous_bitscore_ratio": config.blast_rescue.ambiguous_bitscore_ratio,
                "ambiguous_overlap_ratio": config.blast_rescue.ambiguous_overlap_ratio,
                "its": asdict(config.blast_rescue.its),
                "its2": asdict(config.blast_rescue.its2),
            },
            "tree_shrink_qc": {
                "quantile": config.tree_shrink_qc.quantile,
                "bootstrap": config.tree_shrink_qc.bootstrap,
                "max_removed": (
                    "auto-select"
                    if config.tree_shrink_qc.max_removed is None
                    else config.tree_shrink_qc.max_removed
                ),
            },
        },
    }


def _read_toml(path: Path) -> dict[str, Any]:
    try:
        return tomllib.loads(path.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as error:
        raise ConfigError(f"invalid TOML config {path}: {error}") from error


def _parse_config(raw: dict[str, Any]) -> AppConfig:
    paths = raw.get("paths", {})
    collectors = raw.get("collectors", {})
    genbank = collectors.get("genbank", {})
    build = raw.get("build", {})
    itsxrust = build.get("itsxrust", {})
    blast_rescue = build.get("blast_rescue", {})
    tree_shrink_qc = build.get("tree_shrink_qc", {})

    app_kwargs: dict[str, Any] = {}
    if "data_dir" in paths:
        app_kwargs["data_dir"] = Path(str(paths["data_dir"])).expanduser()

    collector_kwargs = {
        key: collectors[key]
        for key in ("batch_size", "download_workers", "timeout", "retry_attempts")
        if key in collectors
    }
    if "email" in genbank:
        collector_kwargs["genbank_email"] = str(genbank["email"])
    if "api_key" in genbank:
        collector_kwargs["genbank_api_key"] = str(genbank["api_key"] or "") or None

    blast_rescue_kwargs = {
        key: blast_rescue[key]
        for key in (
            "blastn_dust",
            "word_size",
            "evalue",
            "endpoint_margin_bases",
            "endpoint_margin_fraction",
            "ambiguous_bitscore_ratio",
            "ambiguous_overlap_ratio",
        )
        if key in blast_rescue
    }
    if "its" in blast_rescue:
        blast_rescue_kwargs["its"] = replace(BlastRescueConfig().its, **blast_rescue["its"])
    if "its2" in blast_rescue:
        blast_rescue_kwargs["its2"] = replace(BlastRescueConfig().its2, **blast_rescue["its2"])

    itsxrust_kwargs = {
        key: itsxrust[key]
        for key in ("inc_e", "min_anchor_score", "max_per_anchor", "max_anchor_evalue")
        if key in itsxrust
    }

    tree_shrink_kwargs = {
        key: tree_shrink_qc[key]
        for key in ("quantile", "bootstrap")
        if key in tree_shrink_qc
    }
    if "max_removed" in tree_shrink_qc:
        max_removed = tree_shrink_qc["max_removed"]
        tree_shrink_kwargs["max_removed"] = (
            None
            if isinstance(max_removed, str)
            and max_removed.strip().lower() in {"auto", "auto-select"}
            else max_removed
        )

    return AppConfig(
        **app_kwargs,
        collectors=CollectorConfig(**collector_kwargs),
        blast_rescue=BlastRescueConfig(**blast_rescue_kwargs),
        itsxrust=ItsxrustConfig(**itsxrust_kwargs),
        tree_shrink_qc=TreeShrinkConfig(**tree_shrink_kwargs),
    )


def _toml_dump(data: dict[str, Any]) -> str:
    def format_value(value: Any) -> str:
        if isinstance(value, bool):
            return "true" if value else "false"
        if isinstance(value, (int, float)):
            return str(value)
        escaped = str(value).replace("\\", "\\\\").replace('"', '\\"')
        return f'"{escaped}"'

    def write_table(lines: list[str], table: dict[str, Any], prefix: list[str]) -> None:
        scalar_items = [
            (key, value) for key, value in table.items() if not isinstance(value, dict)
        ]
        child_items = [
            (key, value) for key, value in table.items() if isinstance(value, dict)
        ]
        if prefix:
            if lines:
                lines.append("")
            lines.append(f"[{'.'.join(prefix)}]")
        for key, value in scalar_items:
            lines.append(f"{key} = {format_value(value)}")
        for key, value in child_items:
            write_table(lines, value, [*prefix, key])

    lines: list[str] = []
    write_table(lines, data, [])
    return "\n".join(lines).rstrip() + "\n"
