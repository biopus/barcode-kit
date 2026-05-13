from __future__ import annotations

import json
import sys
import threading
import time
from collections import deque
from collections.abc import Iterable
from dataclasses import asdict
from pathlib import Path
from typing import Annotated

import typer

from barcode_kit import config as config_module
from barcode_kit.builder import build_dataset
from barcode_kit.exceptions import BarcodeKitError
from barcode_kit.genbank import SyncService
from barcode_kit.models import ItsExtractionMode, Marker, TaxonQuery
from barcode_kit.phylogeny import TreeShrinkQcConfig
from barcode_kit.storage import Storage
from barcode_kit.taxonomy import ETETaxonomyResolver


app = typer.Typer(help="Build local GenBank-backed DNA barcode datasets.")
db_app = typer.Typer(help="Inspect the local cache database.")
config_app = typer.Typer(help="Inspect or update barcode-kit configuration.")
app.add_typer(db_app, name="db")
app.add_typer(config_app, name="config")
DOWNLOAD_PROGRESS_INTERVAL_SECONDS = 0.5


FamilyOption = Annotated[str | None, typer.Option("--family", help="Family taxon name.")]
GenusOption = Annotated[str | None, typer.Option("--genus", help="Genus taxon name.")]
SpeciesOption = Annotated[str | None, typer.Option("--species", help="Species scientific name.")]
YesOption = Annotated[bool, typer.Option("--yes", "-y", help="Skip confirmation prompt.")]
RankOption = Annotated[
    str | None,
    typer.Option("--rank", help="List cached taxa at rank: family, genus, or species."),
]


@app.command()
def sync(
    marker: Annotated[Marker, typer.Option("--marker", case_sensitive=False)],
    family: FamilyOption = None,
    genus: GenusOption = None,
    species: SpeciesOption = None,
) -> None:
    """Synchronize matching GenBank records into the local cache."""
    _run_user_command(lambda: _sync(marker, family, genus, species))


@app.command()
def build(
    marker: Annotated[Marker, typer.Option("--marker", case_sensitive=False)],
    outdir: Annotated[Path, typer.Option("--outdir")] = Path("."),
    family: FamilyOption = None,
    genus: GenusOption = None,
    species: SpeciesOption = None,
    min_length: Annotated[int | None, typer.Option("--min-length")] = None,
    max_ambiguous_content: Annotated[
        float | None,
        typer.Option("--max-ambiguous-content"),
    ] = None,
    exclude_hybrid: Annotated[bool, typer.Option("--exclude-hybrid")] = False,
    exclude_uncertain: Annotated[bool, typer.Option("--exclude-uncertain")] = False,
    its_extraction_mode: Annotated[
        ItsExtractionMode,
        typer.Option("--its-extraction-mode", case_sensitive=False),
    ] = ItsExtractionMode.HMM_BLAST,
    tree_shrink_qc: Annotated[
        bool,
        typer.Option("--tree-shrink-qc", help="Run MAFFT, IQ-TREE, and TreeShrink long-branch QC."),
    ] = False,
    tree_shrink_qc_threads: Annotated[
        int,
        typer.Option("--tree-shrink-qc-threads", help="Threads for MAFFT and IQ-TREE in TreeShrink QC."),
    ] = 1,
    tree_shrink_qc_quantile: Annotated[
        float,
        typer.Option("--tree-shrink-qc-quantile", help="TreeShrink false positive tolerance quantile."),
    ] = 0.05,
) -> None:
    """Build a FASTA dataset from the local cache."""
    _run_user_command(
        lambda: _build(
            marker,
            outdir,
            family,
            genus,
            species,
            min_length,
            max_ambiguous_content,
            exclude_hybrid,
            exclude_uncertain,
            its_extraction_mode,
            tree_shrink_qc,
            tree_shrink_qc_threads,
            tree_shrink_qc_quantile,
        )
    )


@db_app.command("status")
def db_status() -> None:
    """Show database and cache counts."""
    _run_user_command(lambda: _db_status())


@db_app.command("info")
def db_info(
    rank: RankOption = None,
    family: FamilyOption = None,
    genus: GenusOption = None,
    species: SpeciesOption = None,
) -> None:
    """Show marker coverage or cached taxon summaries."""
    _run_user_command(lambda: _db_info(rank, family, genus, species))


@db_app.command("remove")
def db_remove(
    accession: Annotated[
        str | None,
        typer.Option("--accession", help="Accession root or version to remove."),
    ] = None,
    family: FamilyOption = None,
    genus: GenusOption = None,
    species: SpeciesOption = None,
    yes: YesOption = False,
) -> None:
    """Remove matching records from the local cache."""
    _run_user_command(lambda: _db_remove(accession, family, genus, species, yes))


@db_app.command("clear")
def db_clear(yes: YesOption = False) -> None:
    """Remove all local cache records and GenBank cache files."""
    _run_user_command(lambda: _db_clear(yes))


@db_app.command("prune")
def db_prune(yes: YesOption = False) -> None:
    """Remove inconsistent local cache database rows and files."""
    _run_user_command(lambda: _db_prune(yes))


@config_app.command("list")
def config_list() -> None:
    """Print the effective configuration as JSON."""
    _run_user_command(lambda: _config_list())


@config_app.command("set")
def config_set(key: str, value: str) -> None:
    """Set a configuration value."""
    _run_user_command(lambda: _config_set(key, value))


def main() -> None:
    app()


def _sync(marker: Marker, family: str | None, genus: str | None, species: str | None) -> None:
    query = _taxon_query(family, genus, species)
    config = config_module.load_or_create_config()
    config_module.ensure_app_dirs(config)
    storage = Storage(config.database_path)
    service = SyncService(config, storage, ETETaxonomyResolver())
    progress = DownloadProgressState()
    reporter = TerminalDownloadReporter(
        progress,
        interval_seconds=DOWNLOAD_PROGRESS_INTERVAL_SECONDS,
    )
    try:
        reporter.start()
        result = service.sync(query, marker, progress=progress)
    finally:
        reporter.stop()
        service.close()
    typer.echo(
        json.dumps(asdict(result), ensure_ascii=False, indent=2)
    )
    if result.failed:
        raise typer.Exit(code=2)


class DownloadProgressState:
    def __init__(self, speed_window_seconds: float = 5.0):
        self.speed_window_seconds = speed_window_seconds
        self.started_at: float | None = None
        self.bytes_written = 0
        self.completed_records = 0
        self.failed_records = 0
        self.total_records = 0
        self._lock = threading.Lock()
        self._byte_events: deque[tuple[float, int]] = deque()
        self._recent_bytes = 0

    def start_download(self, total_records: int) -> None:
        with self._lock:
            if self.started_at is None:
                self.started_at = time.monotonic()
            self.total_records = total_records

    def record_downloaded_bytes(self, bytes_delta: int) -> None:
        if bytes_delta <= 0:
            return
        with self._lock:
            self.bytes_written += bytes_delta
            self._byte_events.append((time.monotonic(), bytes_delta))
            self._recent_bytes += bytes_delta
            self._prune_byte_events()

    def record_downloaded_record(self, accession: str) -> None:
        del accession
        with self._lock:
            self.completed_records += 1

    def finish_download(self, succeeded_records: int, failed_records: int) -> None:
        with self._lock:
            self.completed_records = max(self.completed_records, succeeded_records)
            self.failed_records = failed_records

    def snapshot(self) -> dict[str, float | int] | None:
        with self._lock:
            if self.started_at is None:
                return None
            elapsed = max(time.monotonic() - self.started_at, 1e-9)
            self._prune_byte_events()
            return {
                "elapsed": elapsed,
                "bytes_written": self.bytes_written,
                "completed_records": self.completed_records,
                "failed_records": self.failed_records,
                "total_records": self.total_records,
                "average_speed_kb_s": self.bytes_written / elapsed / 1024,
                "recent_speed_kb_s": self._recent_bytes / self.speed_window_seconds / 1024,
            }

    def _prune_byte_events(self) -> None:
        cutoff = time.monotonic() - self.speed_window_seconds
        while self._byte_events and self._byte_events[0][0] < cutoff:
            _, size = self._byte_events.popleft()
            self._recent_bytes -= size


class TerminalDownloadReporter:
    def __init__(self, progress: DownloadProgressState, interval_seconds: float) -> None:
        self.progress = progress
        self.interval_seconds = interval_seconds
        self._stop_event = threading.Event()
        self._thread = threading.Thread(
            target=self._run,
            name="barcode-kit-download-progress",
            daemon=True,
        )

    def start(self) -> None:
        self._thread.start()

    def stop(self) -> None:
        self._stop_event.set()
        self._thread.join()
        self._print_snapshot(final=True)

    def _run(self) -> None:
        while not self._stop_event.wait(self.interval_seconds):
            self._print_snapshot(final=False)

    def _print_snapshot(self, *, final: bool) -> None:
        snapshot = self.progress.snapshot()
        if snapshot is None:
            return
        total_records = int(snapshot["total_records"])
        completed_records = int(snapshot["completed_records"])
        downloaded = (
            f"{completed_records}/{total_records}"
            if total_records > 0
            else str(completed_records)
        )
        line = (
            f"downloaded={downloaded} "
            f"speed={snapshot['recent_speed_kb_s']:.1f} KB/s "
            f"avg={snapshot['average_speed_kb_s']:.1f} KB/s "
            f"bytes={snapshot['bytes_written']} "
            f"failed={snapshot['failed_records']}"
        )
        suffix = "\n" if final else ""
        sys.stderr.write(f"\r\033[2K{line}{suffix}")
        sys.stderr.flush()


def _build(
    marker: Marker,
    outdir: Path,
    family: str | None,
    genus: str | None,
    species: str | None,
    min_length: int | None,
    max_ambiguous_content: float | None,
    exclude_hybrid: bool,
    exclude_uncertain: bool,
    its_extraction_mode: ItsExtractionMode,
    tree_shrink_qc: bool,
    tree_shrink_qc_threads: int,
    tree_shrink_qc_quantile: float,
) -> None:
    query = _taxon_query(family, genus, species)
    config = config_module.load_or_create_config()
    storage = Storage(config.database_path)
    report = build_dataset(
        config,
        storage,
        query,
        marker,
        outdir,
        min_length=min_length,
        max_ambiguous_content=max_ambiguous_content,
        exclude_hybrid=exclude_hybrid,
        exclude_uncertain=exclude_uncertain,
        its_extraction_mode=its_extraction_mode,
        tree_shrink_qc=(
            TreeShrinkQcConfig(
                threads=tree_shrink_qc_threads,
                quantile=tree_shrink_qc_quantile,
            )
            if tree_shrink_qc
            else None
        ),
    )
    included = sum(1 for entry in report if entry.included)
    typer.echo(
        json.dumps(
            {
                "outdir": str(outdir),
                "marker": marker.value,
                "records": len(report),
                "included": included,
                "excluded": len(report) - included,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


def _db_status() -> None:
    config = config_module.load_config()
    storage = Storage(config.database_path)
    storage.initialize()
    payload = {
        "database": str(config.database_path),
        "genbank_cache_dir": str(config.genbank_cache_dir),
        **storage.counts(),
    }
    typer.echo(json.dumps(payload, ensure_ascii=False, indent=2))


def _db_info(
    rank: str | None,
    family: str | None,
    genus: str | None,
    species: str | None,
) -> None:
    config = config_module.load_config()
    storage = Storage(config.database_path)
    storage.initialize()
    query = _taxon_query(family, genus, species) if any([family, genus, species]) else None
    if rank is not None:
        rank = rank.lower()
        _validate_info_rank_filter(rank, query)
        payload = {
            "rank": rank,
            "query": query.__dict__ if query else None,
            "taxa": storage.taxon_summaries(rank, query),
        }
    else:
        payload = {
            "query": query.__dict__ if query else None,
            "markers": storage.marker_counts(query),
        }
    typer.echo(json.dumps(payload, ensure_ascii=False, indent=2))


def _db_remove(
    accession: str | None,
    family: str | None,
    genus: str | None,
    species: str | None,
    yes: bool,
) -> None:
    taxon_selected = any([family, genus, species])
    if bool(accession) == taxon_selected:
        raise BarcodeKitError(
            "provide exactly one of --accession, --family, --genus, or --species"
        )

    config = config_module.load_config()
    storage = Storage(config.database_path)
    storage.initialize()
    query = _taxon_query(family, genus, species) if taxon_selected else None
    records = storage.cache_records(query, accession=accession)
    _confirm_cache_mutation(f"Remove {len(records)} local cache record(s)", yes)

    files_removed = _remove_genbank_cache_files(
        _genbank_cache_path(config.genbank_cache_dir, record.accession_version)
        for record in records
    )
    database_removed = storage.delete_cache_records(record.accession_root for record in records)
    taxonomy_removed = storage.delete_orphan_taxonomy()
    payload = {
        "selector": {"accession": accession, "query": query.__dict__ if query else None},
        "database_records_removed": database_removed,
        "files_removed": files_removed,
        "taxonomy_removed": taxonomy_removed,
    }
    typer.echo(json.dumps(payload, ensure_ascii=False, indent=2))


def _db_clear(yes: bool) -> None:
    config = config_module.load_config()
    storage = Storage(config.database_path)
    storage.initialize()
    records = storage.cache_records()
    cache_files = (
        list(config.genbank_cache_dir.glob("*.gb"))
        if config.genbank_cache_dir.exists()
        else []
    )
    _confirm_cache_mutation(
        f"Remove all {len(records)} local cache record(s) and {len(cache_files)} GenBank file(s)",
        yes,
    )

    files_removed = _remove_genbank_cache_files(cache_files)
    deleted = storage.clear_cache()
    payload = {
        "database_records_removed": deleted["genbank_cache"],
        "files_removed": files_removed,
        "taxonomy_removed": deleted["taxonomy"],
    }
    typer.echo(json.dumps(payload, ensure_ascii=False, indent=2))


def _db_prune(yes: bool) -> None:
    config = config_module.load_config()
    storage = Storage(config.database_path)
    storage.initialize()
    records = storage.cache_records()
    missing_file_records = [
        record
        for record in records
        if not _genbank_cache_path(config.genbank_cache_dir, record.accession_version).exists()
    ]
    cache_files = (
        list(config.genbank_cache_dir.glob("*.gb"))
        if config.genbank_cache_dir.exists()
        else []
    )
    referenced_names = {f"{record.accession_version}.gb" for record in records}
    unreferenced_files = [path for path in cache_files if path.name not in referenced_names]
    _confirm_cache_mutation(
        (
            f"Prune {len(missing_file_records)} database record(s) and "
            f"{len(unreferenced_files)} unreferenced GenBank file(s)"
        ),
        yes,
    )

    files_removed = _remove_genbank_cache_files(unreferenced_files)
    database_removed = storage.delete_cache_records(
        record.accession_root for record in missing_file_records
    )
    taxonomy_removed = storage.delete_orphan_taxonomy()
    payload = {
        "database_records_removed": database_removed,
        "files_removed": files_removed,
        "taxonomy_removed": taxonomy_removed,
    }
    typer.echo(json.dumps(payload, ensure_ascii=False, indent=2))


def _config_list() -> None:
    path = config_module.ensure_config_file()
    config = config_module.load_config(path)
    typer.echo(json.dumps(config_module.config_as_dict(config), ensure_ascii=False, indent=2))


def _config_set(key: str, value: str) -> None:
    config = config_module.set_config_value(key, value)
    typer.echo(json.dumps(config_module.config_as_dict(config), ensure_ascii=False, indent=2))


def _taxon_query(family: str | None, genus: str | None, species: str | None) -> TaxonQuery:
    values = [("family", family), ("genus", genus), ("species", species)]
    selected = [(rank, value) for rank, value in values if value]
    if len(selected) != 1:
        raise BarcodeKitError("provide exactly one of --family, --genus, or --species")
    rank, value = selected[0]
    return TaxonQuery(rank=rank, name=value)


def _validate_info_rank_filter(rank: str, query: TaxonQuery | None) -> None:
    allowed_filters = {
        "family": set(),
        "genus": {"family"},
        "species": {"family", "genus"},
    }
    if rank not in allowed_filters:
        raise BarcodeKitError("rank must be one of family, genus, or species")
    if query is None:
        return
    allowed = allowed_filters[rank]
    if query.rank in allowed:
        return
    if not allowed:
        raise BarcodeKitError(f"--rank {rank} cannot be combined with taxon filters")
    allowed_options = " or ".join(f"--{value}" for value in sorted(allowed))
    raise BarcodeKitError(f"--rank {rank} can only be filtered by {allowed_options}")


def _genbank_cache_path(cache_dir: Path, accession_version: str) -> Path:
    return cache_dir / f"{accession_version}.gb"


def _remove_genbank_cache_files(paths: Iterable[Path]) -> int:
    removed = 0
    for path in dict.fromkeys(paths):
        if not path.exists():
            continue
        path.unlink()
        removed += 1
    return removed


def _confirm_cache_mutation(message: str, yes: bool) -> None:
    if yes:
        return
    if not typer.confirm(f"{message}. Continue?"):
        raise typer.Abort()


def _run_user_command(callback) -> None:
    try:
        callback()
    except BarcodeKitError as error:
        typer.echo(f"Error: {error}", err=True)
        raise typer.Exit(code=1) from error


if __name__ == "__main__":
    main()
