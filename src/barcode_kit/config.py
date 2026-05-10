from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from barcode_kit.exceptions import ConfigError


DEFAULT_DATA_DIR = Path.home() / ".barcode-kit"
DEFAULT_CONFIG = {
    "paths": {"data_dir": str(DEFAULT_DATA_DIR)},
    "collectors": {
        "batch_size": 500,
        "download_workers": 8,
        "timeout": 30,
        "retry_attempts": 3,
        "genbank": {
            "email": "",
            "api_key": "",
        },
    },
    "build": {
        "itsxrust": {
            "inc_e": 0.01,
            "min_anchor_score": 8,
            "max_per_anchor": 20,
            "max_anchor_evalue": 0.01,
        },
        "blast_rescue": {
            "blastn_dust": "no",
            "word_size": 11,
            "evalue": 1e-3,
            "endpoint_margin_bases": 15,
            "endpoint_margin_fraction": 0.05,
            "ambiguous_bitscore_ratio": 0.95,
            "ambiguous_overlap_ratio": 0.80,
            "its": {
                "min_subject_coverage": 0.85,
                "min_identity": 0.80,
                "min_query_length_ratio": 0.85,
                "max_query_length_ratio": 1.20,
            },
            "its2": {
                "min_subject_coverage": 0.90,
                "min_identity": 0.85,
                "min_query_length_ratio": 0.90,
                "max_query_length_ratio": 1.15,
            },
        },
    },
}


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
    endpoint_margin_bases: int = 30
    endpoint_margin_fraction: float = 0.05
    ambiguous_bitscore_ratio: float = 0.95
    ambiguous_overlap_ratio: float = 0.80
    its: BlastRescueMarkerConfig = field(
        default_factory=lambda: BlastRescueMarkerConfig(
            min_subject_coverage=0.85,
            min_identity=0.80,
            min_query_length_ratio=0.75,
            max_query_length_ratio=1.30,
        )
    )
    its2: BlastRescueMarkerConfig = field(
        default_factory=lambda: BlastRescueMarkerConfig(
            min_subject_coverage=0.90,
            min_identity=0.85,
            min_query_length_ratio=0.70,
            max_query_length_ratio=1.35,
        )
    )


@dataclass(frozen=True)
class ItsxrustConfig:
    inc_e: float = 0.01
    min_anchor_score: int = 8
    max_per_anchor: int = 20
    max_anchor_evalue: float = 0.01


@dataclass(frozen=True)
class AppConfig:
    data_dir: Path
    batch_size: int
    download_workers: int
    timeout: float
    retry_attempts: int
    genbank_email: str
    genbank_api_key: str | None = None
    blast_rescue: BlastRescueConfig = field(default_factory=BlastRescueConfig)
    itsxrust: ItsxrustConfig = field(default_factory=ItsxrustConfig)

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


def default_config_path() -> Path:
    return Path(os.environ.get("BARCODE_KIT_CONFIG", DEFAULT_DATA_DIR / "config.toml"))


def load_config(path: Path | None = None) -> AppConfig:
    path = path or default_config_path()
    if path.exists():
        raw = _deep_merge(DEFAULT_CONFIG, _read_toml(path))
    else:
        raw = DEFAULT_CONFIG
    return _parse_config(raw)


def load_or_create_config(path: Path | None = None) -> AppConfig:
    path = ensure_config_file(path)
    return load_config(path)


def ensure_app_dirs(config: AppConfig) -> None:
    config.data_dir.mkdir(parents=True, exist_ok=True)
    config.genbank_cache_dir.mkdir(parents=True, exist_ok=True)
    config.logs_dir.mkdir(parents=True, exist_ok=True)


def ensure_config_file(path: Path | None = None) -> Path:
    path = path or default_config_path()
    config = load_config(path)
    ensure_app_dirs(config)
    if not path.exists():
        write_config(config, path)
    return path


def write_config(config: AppConfig, path: Path | None = None) -> None:
    path = path or config.config_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_toml_dump(config_as_dict(config)), encoding="utf-8")


def set_config_value(key: str, value: str, path: Path | None = None) -> AppConfig:
    path = path or default_config_path()
    raw = _deep_merge(DEFAULT_CONFIG, _read_toml(path) if path.exists() else {})
    parts = key.split(".")
    parts_tuple = tuple(parts)
    if parts_tuple not in {
        ("paths", "data_dir"),
        ("collectors", "batch_size"),
        ("collectors", "download_workers"),
        ("collectors", "timeout"),
        ("collectors", "retry_attempts"),
        ("collectors", "genbank", "email"),
        ("collectors", "genbank", "api_key"),
        ("genbank", "email"),
        ("genbank", "api_key"),
        ("build", "itsxrust", "inc_e"),
        ("build", "itsxrust", "min_anchor_score"),
        ("build", "itsxrust", "max_per_anchor"),
        ("build", "itsxrust", "max_anchor_evalue"),
        ("build", "blast_rescue", "blastn_dust"),
        ("build", "blast_rescue", "word_size"),
        ("build", "blast_rescue", "evalue"),
        ("build", "blast_rescue", "endpoint_margin_bases"),
        ("build", "blast_rescue", "endpoint_margin_fraction"),
        ("build", "blast_rescue", "ambiguous_bitscore_ratio"),
        ("build", "blast_rescue", "ambiguous_overlap_ratio"),
        ("build", "blast_rescue", "its", "min_subject_coverage"),
        ("build", "blast_rescue", "its", "min_identity"),
        ("build", "blast_rescue", "its", "min_query_length_ratio"),
        ("build", "blast_rescue", "its", "max_query_length_ratio"),
        ("build", "blast_rescue", "its2", "min_subject_coverage"),
        ("build", "blast_rescue", "its2", "min_identity"),
        ("build", "blast_rescue", "its2", "min_query_length_ratio"),
        ("build", "blast_rescue", "its2", "max_query_length_ratio"),
    }:
        raise ConfigError(f"unknown config key: {key}")
    if parts[0] == "genbank":
        parts = ["collectors", "genbank", *parts[1:]]
    current: dict[str, Any] = raw
    for part in parts[:-1]:
        current = current.setdefault(part, {})
    current[parts[-1]] = _coerce_config_value(key, value)
    config = _parse_config(raw)
    ensure_app_dirs(config)
    write_config(config, path)
    return config


def config_as_dict(config: AppConfig) -> dict[str, Any]:
    return {
        "paths": {"data_dir": str(config.data_dir)},
        "collectors": {
            "batch_size": config.batch_size,
            "download_workers": config.download_workers,
            "timeout": config.timeout,
            "retry_attempts": config.retry_attempts,
            "genbank": {
                "email": config.genbank_email,
                "api_key": config.genbank_api_key or "",
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
                "its": _blast_marker_config_as_dict(config.blast_rescue.its),
                "its2": _blast_marker_config_as_dict(config.blast_rescue.its2),
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
    return AppConfig(
        data_dir=Path(str(paths.get("data_dir") or DEFAULT_DATA_DIR)).expanduser(),
        batch_size=int(collectors.get("batch_size", 500)),
        download_workers=int(collectors.get("download_workers", 1)),
        timeout=float(collectors.get("timeout", 30)),
        retry_attempts=int(collectors.get("retry_attempts", 3)),
        genbank_email=str(genbank.get("email", "")),
        genbank_api_key=str(genbank.get("api_key") or "") or None,
        blast_rescue=_parse_blast_rescue_config(blast_rescue),
        itsxrust=_parse_itsxrust_config(itsxrust),
    )


def _coerce_config_value(key: str, value: str) -> str | int | float:
    if key.endswith(
        (
            "batch_size",
            "download_workers",
            "retry_attempts",
            "endpoint_margin_bases",
            "word_size",
            "min_anchor_score",
            "max_per_anchor",
        )
    ):
        return int(value)
    if key.endswith(
        (
            "timeout",
            "evalue",
            "inc_e",
            "max_anchor_evalue",
            "endpoint_margin_fraction",
            "ambiguous_bitscore_ratio",
            "ambiguous_overlap_ratio",
            "min_subject_coverage",
            "min_identity",
            "min_query_length_ratio",
            "max_query_length_ratio",
        )
    ):
        return float(value)
    return value


def _parse_blast_rescue_config(raw: dict[str, Any]) -> BlastRescueConfig:
    return BlastRescueConfig(
        blastn_dust=str(raw.get("blastn_dust", "no")),
        word_size=int(raw.get("word_size", 11)),
        evalue=float(raw.get("evalue", 1e-3)),
        endpoint_margin_bases=int(raw.get("endpoint_margin_bases", 30)),
        endpoint_margin_fraction=float(raw.get("endpoint_margin_fraction", 0.05)),
        ambiguous_bitscore_ratio=float(raw.get("ambiguous_bitscore_ratio", 0.95)),
        ambiguous_overlap_ratio=float(raw.get("ambiguous_overlap_ratio", 0.80)),
        its=_parse_blast_marker_config(
            raw.get("its", {}),
            BlastRescueMarkerConfig(0.85, 0.80, 0.75, 1.30),
        ),
        its2=_parse_blast_marker_config(
            raw.get("its2", {}),
            BlastRescueMarkerConfig(0.90, 0.85, 0.70, 1.35),
        ),
    )


def _parse_itsxrust_config(raw: dict[str, Any]) -> ItsxrustConfig:
    return ItsxrustConfig(
        inc_e=float(raw.get("inc_e", 0.01)),
        min_anchor_score=int(raw.get("min_anchor_score", 8)),
        max_per_anchor=int(raw.get("max_per_anchor", 20)),
        max_anchor_evalue=float(raw.get("max_anchor_evalue", 0.01)),
    )


def _parse_blast_marker_config(
    raw: dict[str, Any],
    defaults: BlastRescueMarkerConfig,
) -> BlastRescueMarkerConfig:
    return BlastRescueMarkerConfig(
        min_subject_coverage=float(
            raw.get("min_subject_coverage", defaults.min_subject_coverage)
        ),
        min_identity=float(raw.get("min_identity", defaults.min_identity)),
        min_query_length_ratio=float(
            raw.get("min_query_length_ratio", defaults.min_query_length_ratio)
        ),
        max_query_length_ratio=float(
            raw.get("max_query_length_ratio", defaults.max_query_length_ratio)
        ),
    )


def _blast_marker_config_as_dict(config: BlastRescueMarkerConfig) -> dict[str, float]:
    return {
        "min_subject_coverage": config.min_subject_coverage,
        "min_identity": config.min_identity,
        "min_query_length_ratio": config.min_query_length_ratio,
        "max_query_length_ratio": config.max_query_length_ratio,
    }


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in base.items():
        if isinstance(value, dict):
            result[key] = _deep_merge(value, override.get(key, {}))
        else:
            result[key] = value
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def _toml_dump(data: dict[str, Any]) -> str:
    lines: list[str] = []
    _write_table(lines, data, [])
    return "\n".join(lines).rstrip() + "\n"


def _write_table(lines: list[str], table: dict[str, Any], prefix: list[str]) -> None:
    scalar_items = [(key, value) for key, value in table.items() if not isinstance(value, dict)]
    child_items = [(key, value) for key, value in table.items() if isinstance(value, dict)]
    if prefix:
        if lines:
            lines.append("")
        lines.append(f"[{'.'.join(prefix)}]")
    for key, value in scalar_items:
        lines.append(f"{key} = {_format_toml_value(value)}")
    for key, value in child_items:
        _write_table(lines, value, [*prefix, key])


def _format_toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    escaped = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'
