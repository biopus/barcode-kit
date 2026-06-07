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
from barcode_kit import validation as validation_module
from barcode_kit.builder import build_dataset
from barcode_kit.exceptions import BarcodeKitError, GenBankError, TaxonomyError
from barcode_kit.genbank import SyncService
from barcode_kit.models import GenBankCacheRecord, Marker, TaxonConstraint, TaxonExclusion, TaxonQuery
from barcode_kit.parser import parse_genbank_file
from barcode_kit.storage import Storage
from barcode_kit.taxonomy import ETETaxonomyResolver


__all__ = [
    "DownloadProgressState",
    "TerminalDownloadReporter",
    "app",
    "build",
    "db_app",
    "db_clear",
    "db_info",
    "db_prune",
    "db_rebuild",
    "db_remove",
    "db_status",
    "main",
    "qc",
    "sync",
]


app = typer.Typer(help="Build local GenBank-backed DNA barcode datasets.")
db_app = typer.Typer(help="Inspect the local cache database.")
app.add_typer(db_app, name="db")
DOWNLOAD_PROGRESS_INTERVAL_SECONDS = 0.5


KingdomOption = Annotated[
    str | None,
    typer.Option("--kingdom", "--kindom", help="Kingdom taxon name."),
]
PhylumOption = Annotated[str | None, typer.Option("--phylum", help="Phylum taxon name.")]
ClassOption = Annotated[str | None, typer.Option("--class", help="Class taxon name.")]
OrderOption = Annotated[str | None, typer.Option("--order", help="Order taxon name.")]
FamilyOption = Annotated[str | None, typer.Option("--family", help="Family taxon name.")]
GenusOption = Annotated[str | None, typer.Option("--genus", help="Genus taxon name.")]
SpeciesOption = Annotated[str | None, typer.Option("--species", help="Species scientific name.")]
TaxidOption = Annotated[int | None, typer.Option("--taxid", help="NCBI taxonomy ID.")]
YesOption = Annotated[bool, typer.Option("--yes", "-y", help="Skip confirmation prompt.")]
RankOption = Annotated[
    str | None,
    typer.Option("--rank", help="List cached taxa at rank: family, genus, or species."),
]
ExcludeOption = Annotated[
    list[TaxonExclusion] | None,
    typer.Option("--exclude", case_sensitive=False, help="Exclude taxon class from QC output."),
]


@app.command()
def sync(
    marker: Annotated[Marker, typer.Option("--marker", case_sensitive=False)],
    taxid: TaxidOption = None,
    kingdom: KingdomOption = None,
    phylum: PhylumOption = None,
    class_name: ClassOption = None,
    order: OrderOption = None,
    family: FamilyOption = None,
    genus: GenusOption = None,
    species: SpeciesOption = None,
) -> None:
    """Synchronize matching GenBank records into the local cache."""
    _run_user_command(
        lambda: _sync(
            marker,
            taxid,
            kingdom,
            phylum,
            class_name,
            order,
            family,
            genus,
            species,
        )
    )


@app.command()
def build(
    marker: Annotated[Marker, typer.Option("--marker", case_sensitive=False)],
    outdir: Annotated[Path, typer.Option("--outdir")] = Path("."),
    kingdom: KingdomOption = None,
    phylum: PhylumOption = None,
    class_name: ClassOption = None,
    order: OrderOption = None,
    family: FamilyOption = None,
    genus: GenusOption = None,
    species: SpeciesOption = None,
) -> None:
    """Build a raw FASTA dataset from the local cache."""
    _run_user_command(
        lambda: _build(
            marker,
            outdir,
            kingdom,
            phylum,
            class_name,
            order,
            family,
            genus,
            species,
        )
    )


@app.command()
def qc(
    dataset: Annotated[Path, typer.Option("--dataset")],
    min_length: Annotated[int | None, typer.Option("--min-length")] = None,
    max_ambiguous_content: Annotated[
        float | None,
        typer.Option("--max-ambiguous"),
    ] = None,
    exclude: ExcludeOption = None,
    enable_tree_shrink_qc: Annotated[
        bool,
        typer.Option("--tree-shrink-qc", help="Run MAFFT, IQ-TREE, and TreeShrink long-branch QC."),
    ] = False,
) -> None:
    """Run selected QC checks on a barcode-kit dataset."""
    _run_user_command(
        lambda: _qc(
            dataset,
            min_length,
            max_ambiguous_content,
            exclude,
            enable_tree_shrink_qc,
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
    marker: Annotated[
        Marker | None,
        typer.Option("--marker", case_sensitive=False, help="Marker records to remove."),
    ] = None,
    kingdom: KingdomOption = None,
    phylum: PhylumOption = None,
    class_name: ClassOption = None,
    order: OrderOption = None,
    family: FamilyOption = None,
    genus: GenusOption = None,
    species: SpeciesOption = None,
    yes: YesOption = False,
) -> None:
    """Remove matching records from the local cache."""
    _run_user_command(
        lambda: _db_remove(
            accession,
            marker,
            kingdom,
            phylum,
            class_name,
            order,
            family,
            genus,
            species,
            yes,
        )
    )


@db_app.command("clear")
def db_clear(
    yes: YesOption = False,
    metadata_only: Annotated[
        bool,
        typer.Option("--metadata-only", help="Keep GenBank cache files and remove only metadata."),
    ] = False,
) -> None:
    """Remove all local cache records and GenBank cache files."""
    _run_user_command(lambda: _db_clear(yes, metadata_only))


@db_app.command("rebuild")
def db_rebuild(yes: YesOption = False) -> None:
    """Rebuild local metadata from GenBank cache files."""
    _run_user_command(lambda: _db_rebuild(yes))


@db_app.command("prune")
def db_prune(yes: YesOption = False) -> None:
    """Remove inconsistent local cache database rows and files."""
    _run_user_command(lambda: _db_prune(yes))


def main() -> None:
    app()


def _sync(
    marker: Marker,
    taxid: int | None,
    kingdom: str | None,
    phylum: str | None,
    class_name: str | None,
    order: str | None,
    family: str | None,
    genus: str | None,
    species: str | None,
) -> None:
    query = _sync_query(taxid, kingdom, phylum, class_name, order, family, genus, species)
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
    _echo_json(asdict(result))
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

    def record_downloaded_record(self, _accession: str) -> None:
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
    kingdom: str | None,
    phylum: str | None,
    class_name: str | None,
    order: str | None,
    family: str | None,
    genus: str | None,
    species: str | None,
) -> None:
    query = _constrained_taxon_query(
        kingdom,
        phylum,
        class_name,
        order,
        family,
        genus,
        species,
    )
    config = config_module.load_or_create_config()
    storage = Storage(config.database_path)
    report = build_dataset(
        config,
        storage,
        query,
        marker,
        outdir,
    )
    included = sum(1 for entry in report if entry.included)
    _echo_json(
        {
            "outdir": str(outdir),
            "marker": marker.value,
            "records": len(report),
            "included": included,
            "excluded": len(report) - included,
        }
    )


def _qc(
    dataset: Path,
    min_length: int | None,
    max_ambiguous_content: float | None,
    exclude: list[TaxonExclusion] | None,
    enable_tree_shrink_qc: bool,
) -> None:
    selected_exclusions = set(exclude or [])
    if (
        not selected_exclusions
        and min_length is None
        and max_ambiguous_content is None
        and not enable_tree_shrink_qc
    ):
        raise BarcodeKitError("select at least one QC option")

    config = config_module.load_or_create_config()
    qc_report = validation_module.run_qc(
        dataset,
        exclude=selected_exclusions,
        min_length=min_length,
        max_ambiguous_content=max_ambiguous_content,
        enable_tree_shrink_qc=enable_tree_shrink_qc,
        tree_shrink_config=config.tree_shrink_qc,
    )
    records = qc_report["records"]
    included = sum(1 for record in records if record["included"])
    _echo_json(
        {
            "dataset": str(dataset),
            "records": len(records),
            "included": included,
            "excluded": len(records) - included,
        }
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
    _echo_json(payload)


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
            "query": _query_payload(query),
            "taxa": storage.taxon_summaries(rank, query),
        }
    else:
        payload = {
            "query": _query_payload(query),
            "markers": storage.marker_counts(query),
        }
    _echo_json(payload)


def _db_remove(
    accession: str | None,
    marker: Marker | None,
    kingdom: str | None,
    phylum: str | None,
    class_name: str | None,
    order: str | None,
    family: str | None,
    genus: str | None,
    species: str | None,
    yes: bool,
) -> None:
    taxon_selected = any([kingdom, phylum, class_name, order, family, genus, species])
    if not any([accession, marker, taxon_selected]):
        raise BarcodeKitError("provide --accession, --marker, or a taxon filter")

    config = config_module.load_config()
    storage = Storage(config.database_path)
    storage.initialize()
    query = (
        _constrained_taxon_query(
            kingdom,
            phylum,
            class_name,
            order,
            family,
            genus,
            species,
        )
        if taxon_selected
        else None
    )
    records = storage.cache_records(query, accession=accession, marker=marker)
    _confirm_cache_mutation(f"Remove {len(records)} local cache record(s)", yes)

    files_removed = _remove_genbank_cache_files(
        _genbank_cache_path(config.genbank_cache_dir, record.accession_version)
        for record in records
    )
    database_removed = storage.delete_cache_records(record.accession_root for record in records)
    taxonomy_removed = storage.delete_orphan_taxonomy()
    payload = {
        "selector": {
            "accession": accession,
            "marker": marker.value if marker is not None else None,
            "query": _query_payload(query),
        },
        "database_records_removed": database_removed,
        "files_removed": files_removed,
        "taxonomy_removed": taxonomy_removed,
    }
    _echo_json(payload)


def _db_clear(yes: bool, metadata_only: bool = False) -> None:
    config = config_module.load_config()
    storage = Storage(config.database_path)
    storage.initialize()
    records = storage.cache_records()
    cache_files = (
        list(config.genbank_cache_dir.glob("*.gb"))
        if config.genbank_cache_dir.exists()
        else []
    )
    message = f"Remove all {len(records)} local cache record(s)"
    if not metadata_only:
        message += f" and {len(cache_files)} GenBank file(s)"
    _confirm_cache_mutation(message, yes)

    files_removed = 0 if metadata_only else _remove_genbank_cache_files(cache_files)
    deleted = storage.clear_cache()
    payload = {
        "database_records_removed": deleted["genbank_cache"],
        "files_removed": files_removed,
        "taxonomy_removed": deleted["taxonomy"],
    }
    _echo_json(payload)


def _db_rebuild(yes: bool) -> None:
    config = config_module.load_config()
    storage = Storage(config.database_path)
    storage.initialize()
    cache_files = (
        sorted(config.genbank_cache_dir.glob("*.gb"))
        if config.genbank_cache_dir.exists()
        else []
    )
    _confirm_cache_mutation(f"Rebuild metadata from {len(cache_files)} GenBank file(s)", yes)

    resolver = ETETaxonomyResolver()
    rebuilt = 0
    failed: list[dict[str, str]] = []
    for path in cache_files:
        try:
            parsed = parse_genbank_file(path)
            taxonomy = resolver.standardize(parsed.organism, parsed.taxon_id)
            cache_record = GenBankCacheRecord(
                accession_root=parsed.accession.root,
                version=parsed.accession.version,
                accession_version=parsed.accession.value,
                taxon_id=taxonomy.taxon_id,
                has_its=parsed.marker_flags[Marker.ITS],
                has_matk=parsed.marker_flags[Marker.MATK],
                has_rbcl=parsed.marker_flags[Marker.RBCL],
                has_its2=parsed.marker_flags[Marker.ITS2],
            )
            with storage.connect() as connection:
                storage.upsert_taxonomy(taxonomy, connection)
                storage.upsert_genbank_cache(cache_record, connection)
            rebuilt += 1
        except (GenBankError, TaxonomyError) as error:
            failed.append({"file": str(path), "error": str(error)})

    _echo_json(
        {
            "files_scanned": len(cache_files),
            "records_rebuilt": rebuilt,
            "failed": failed,
        }
    )


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
    _echo_json(payload)


def _echo_json(payload: object) -> None:
    typer.echo(json.dumps(payload, ensure_ascii=False, indent=2))


LINEAGE_OPTION_ORDER = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
PRIMARY_TAXON_RANKS = ("family", "genus", "species")


def _sync_query(
    taxid: int | None,
    kingdom: str | None,
    phylum: str | None,
    class_name: str | None,
    order: str | None,
    family: str | None,
    genus: str | None,
    species: str | None,
) -> TaxonQuery:
    constraints = _lineage_constraints(
        kingdom,
        phylum,
        class_name,
        order,
        family,
        genus,
        species,
    )
    if taxid is not None:
        return TaxonQuery(rank="taxid", name=str(taxid), constraints=constraints)
    if constraints:
        return _constrained_taxon_query(
            kingdom,
            phylum,
            class_name,
            order,
            family,
            genus,
            species,
        )
    raise BarcodeKitError("provide --taxid or one of --family, --genus, or --species")


def _constrained_taxon_query(
    kingdom: str | None,
    phylum: str | None,
    class_name: str | None,
    order: str | None,
    family: str | None,
    genus: str | None,
    species: str | None,
) -> TaxonQuery:
    selected = _selected_lineage_values(
        kingdom,
        phylum,
        class_name,
        order,
        family,
        genus,
        species,
    )
    primary = next(
        ((rank, name) for rank, name in reversed(selected) if rank in PRIMARY_TAXON_RANKS),
        None,
    )
    if primary is None:
        raise BarcodeKitError("provide one of --family, --genus, or --species as the target taxon")
    primary_rank, primary_name = primary
    constraints = tuple(
        TaxonConstraint(rank=rank, name=name)
        for rank, name in selected
        if rank != primary_rank
    )
    return TaxonQuery(rank=primary_rank, name=primary_name, constraints=constraints)


def _lineage_constraints(
    kingdom: str | None,
    phylum: str | None,
    class_name: str | None,
    order: str | None,
    family: str | None,
    genus: str | None,
    species: str | None,
) -> tuple[TaxonConstraint, ...]:
    return tuple(
        TaxonConstraint(rank=rank, name=name)
        for rank, name in _selected_lineage_values(
            kingdom,
            phylum,
            class_name,
            order,
            family,
            genus,
            species,
        )
    )


def _selected_lineage_values(
    kingdom: str | None,
    phylum: str | None,
    class_name: str | None,
    order: str | None,
    family: str | None,
    genus: str | None,
    species: str | None,
) -> list[tuple[str, str]]:
    raw = {
        "kingdom": kingdom,
        "phylum": phylum,
        "class": class_name,
        "order": order,
        "family": family,
        "genus": genus,
        "species": species,
    }
    return [
        (rank, value.strip())
        for rank in LINEAGE_OPTION_ORDER
        if (value := raw[rank]) is not None and value.strip()
    ]


def _query_payload(query: TaxonQuery | None) -> dict[str, object] | None:
    if query is None:
        return None
    payload: dict[str, object] = {"rank": query.rank, "name": query.name}
    if query.constraints:
        payload["constraints"] = [
            {"rank": constraint.rank, "name": constraint.name}
            for constraint in query.constraints
        ]
    return payload


def _taxon_query(
    family: str | None,
    genus: str | None,
    species: str | None,
) -> TaxonQuery:
    values = [
        ("family", family),
        ("genus", genus),
        ("species", species),
    ]
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
