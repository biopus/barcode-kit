from __future__ import annotations

import sqlite3
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TypeVar

from barcode_kit.models import GenBankCacheRecord, Marker, TaxonQuery, TaxonomyRecord


__all__ = ["Storage"]


SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS taxonomy (
    taxon_id INTEGER PRIMARY KEY,
    scientific_name TEXT NOT NULL,
    kingdom TEXT,
    phylum TEXT,
    class TEXT,
    "order" TEXT,
    family TEXT,
    genus TEXT,
    species TEXT,
    infraspecific_rank TEXT,
    is_hybrid INTEGER NOT NULL DEFAULT 0,
    is_uncertain INTEGER NOT NULL DEFAULT 0
);

CREATE INDEX IF NOT EXISTS idx_taxonomy_scientific_name ON taxonomy(scientific_name);
CREATE INDEX IF NOT EXISTS idx_taxonomy_kingdom ON taxonomy(kingdom);
CREATE INDEX IF NOT EXISTS idx_taxonomy_phylum ON taxonomy(phylum);
CREATE INDEX IF NOT EXISTS idx_taxonomy_class ON taxonomy("class");
CREATE INDEX IF NOT EXISTS idx_taxonomy_order ON taxonomy("order");
CREATE INDEX IF NOT EXISTS idx_taxonomy_family ON taxonomy(family);
CREATE INDEX IF NOT EXISTS idx_taxonomy_genus ON taxonomy(genus);

CREATE TABLE IF NOT EXISTS genbank_cache (
    accession_root TEXT PRIMARY KEY,
    version INTEGER NOT NULL,
    accession_version TEXT NOT NULL UNIQUE,
    taxon_id INTEGER NOT NULL REFERENCES taxonomy(taxon_id),
    has_its INTEGER NOT NULL DEFAULT 0,
    has_matk INTEGER NOT NULL DEFAULT 0,
    has_rbcl INTEGER NOT NULL DEFAULT 0,
    has_its2 INTEGER NOT NULL DEFAULT 0,
    updated_at TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_genbank_cache_accession_version ON genbank_cache(accession_version);
CREATE INDEX IF NOT EXISTS idx_genbank_cache_taxon_id ON genbank_cache(taxon_id);
"""

SQLITE_PARAMETER_CHUNK_SIZE = 900
TAXONOMY_FILTER_COLUMNS = {
    "kingdom": "kingdom",
    "phylum": "phylum",
    "class": '"class"',
    "order": '"order"',
    "family": "family",
    "genus": "genus",
    "species": "scientific_name",
}
T = TypeVar("T")


class Storage:
    def __init__(self, path: Path):
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.path)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        return connection

    def initialize(self) -> None:
        with self.connect() as connection:
            connection.executescript(SCHEMA)

    def upsert_taxonomy(self, record: TaxonomyRecord, connection: sqlite3.Connection | None = None) -> None:
        sql = """
        INSERT INTO taxonomy (
            taxon_id, scientific_name, kingdom, phylum, "class", "order", family, genus,
            species, infraspecific_rank, is_hybrid, is_uncertain
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(taxon_id) DO UPDATE SET
            scientific_name = excluded.scientific_name,
            kingdom = excluded.kingdom,
            phylum = excluded.phylum,
            "class" = excluded."class",
            "order" = excluded."order",
            family = excluded.family,
            genus = excluded.genus,
            species = excluded.species,
            infraspecific_rank = excluded.infraspecific_rank,
            is_hybrid = excluded.is_hybrid,
            is_uncertain = excluded.is_uncertain
        """
        params = (
            record.taxon_id,
            record.scientific_name,
            record.kingdom,
            record.phylum,
            record.class_name,
            record.order,
            record.family,
            record.genus,
            record.species,
            record.infraspecific_rank,
            int(record.is_hybrid),
            int(record.is_uncertain),
        )
        if connection is not None:
            connection.execute(sql, params)
            return
        with self.connect() as owned:
            owned.execute(sql, params)

    def upsert_genbank_cache(
        self,
        record: GenBankCacheRecord,
        connection: sqlite3.Connection | None = None,
    ) -> None:
        updated_at = (record.updated_at or datetime.now(timezone.utc)).isoformat()
        sql = """
        INSERT INTO genbank_cache (
            accession_root, version, accession_version, taxon_id, has_its, has_matk,
            has_rbcl, has_its2, updated_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(accession_root) DO UPDATE SET
            version = excluded.version,
            accession_version = excluded.accession_version,
            taxon_id = excluded.taxon_id,
            has_its = excluded.has_its,
            has_matk = excluded.has_matk,
            has_rbcl = excluded.has_rbcl,
            has_its2 = excluded.has_its2,
            updated_at = excluded.updated_at
        """
        params = (
            record.accession_root,
            record.version,
            record.accession_version,
            record.taxon_id,
            int(record.has_its),
            int(record.has_matk),
            int(record.has_rbcl),
            int(record.has_its2),
            updated_at,
        )
        if connection is not None:
            connection.execute(sql, params)
            return
        with self.connect() as owned:
            owned.execute(sql, params)

    def get_cache_by_roots(self, accession_roots: Iterable[str]) -> dict[str, GenBankCacheRecord]:
        roots = list(accession_roots)
        if not roots:
            return {}
        placeholders = ",".join("?" for _ in roots)
        with self.connect() as connection:
            rows = connection.execute(
                f"SELECT * FROM genbank_cache WHERE accession_root IN ({placeholders})",
                roots,
            ).fetchall()
        return {row["accession_root"]: _cache_from_row(row) for row in rows}

    def get_cached_versions_by_roots(self, accession_roots: Iterable[str]) -> dict[str, int]:
        roots = list(dict.fromkeys(accession_roots))
        if not roots:
            return {}

        versions: dict[str, int] = {}
        with self.connect() as connection:
            for chunk in _chunks(roots, SQLITE_PARAMETER_CHUNK_SIZE):
                placeholders = ",".join("?" for _ in chunk)
                rows = connection.execute(
                    f"""
                    SELECT accession_root, version
                    FROM genbank_cache
                    WHERE accession_root IN ({placeholders})
                    """,
                    chunk,
                ).fetchall()
                versions.update({row["accession_root"]: int(row["version"]) for row in rows})
        return versions

    def cache_records(
        self,
        query: TaxonQuery | None = None,
        *,
        accession: str | None = None,
        marker: Marker | None = None,
    ) -> list[GenBankCacheRecord]:
        joins = ""
        where: list[str] = []
        params: list[Any] = []
        if query is not None:
            joins = "JOIN taxonomy t ON t.taxon_id = c.taxon_id"
            query_where, query_params = _query_filters(query, "t")
            where.extend(query_where)
            params.extend(query_params)
        if accession is not None:
            where.append("(c.accession_root = ? OR c.accession_version = ?)")
            params.extend([accession, accession])
        if marker is not None:
            where.append(f"c.{marker.cache_column} = 1")

        where_sql = f"WHERE {' AND '.join(where)}" if where else ""
        sql = f"""
        SELECT c.*
        FROM genbank_cache c
        {joins}
        {where_sql}
        ORDER BY c.accession_version
        """
        with self.connect() as connection:
            rows = connection.execute(sql, params).fetchall()
        return [_cache_from_row(row) for row in rows]

    def delete_cache_records(self, accession_roots: Iterable[str]) -> int:
        roots = list(dict.fromkeys(accession_roots))
        if not roots:
            return 0

        deleted = 0
        with self.connect() as connection:
            for chunk in _chunks(roots, SQLITE_PARAMETER_CHUNK_SIZE):
                placeholders = ",".join("?" for _ in chunk)
                cursor = connection.execute(
                    f"DELETE FROM genbank_cache WHERE accession_root IN ({placeholders})",
                    chunk,
                )
                deleted += int(cursor.rowcount)
        return deleted

    def delete_orphan_taxonomy(self) -> int:
        with self.connect() as connection:
            cursor = connection.execute(
                """
                DELETE FROM taxonomy
                WHERE NOT EXISTS (
                    SELECT 1 FROM genbank_cache c WHERE c.taxon_id = taxonomy.taxon_id
                )
                """
            )
        return int(cursor.rowcount)

    def clear_cache(self) -> dict[str, int]:
        with self.connect() as connection:
            cache_cursor = connection.execute("DELETE FROM genbank_cache")
            taxonomy_cursor = connection.execute("DELETE FROM taxonomy")
        return {
            "genbank_cache": int(cache_cursor.rowcount),
            "taxonomy": int(taxonomy_cursor.rowcount),
        }

    def candidate_records(self, query: TaxonQuery, marker: Marker) -> list[tuple[GenBankCacheRecord, TaxonomyRecord]]:
        marker_column = marker.cache_column
        query_where, query_params = _query_filters(query, "t")
        sql = f"""
        SELECT
            c.*,
            t.taxon_id AS t_taxon_id,
            t.scientific_name,
            t.kingdom,
            t.phylum,
            t."class",
            t."order",
            t.family,
            t.genus,
            t.species,
            t.infraspecific_rank,
            t.is_hybrid,
            t.is_uncertain
        FROM genbank_cache c
        JOIN taxonomy t ON t.taxon_id = c.taxon_id
        WHERE c.{marker_column} = 1 AND {' AND '.join(query_where)}
        ORDER BY t.scientific_name, c.accession_version
        """
        with self.connect() as connection:
            rows = connection.execute(sql, query_params).fetchall()
        return [(_cache_from_row(row), _taxonomy_from_joined_row(row)) for row in rows]

    def counts(self) -> dict[str, int]:
        with self.connect() as connection:
            taxonomy = connection.execute("SELECT COUNT(*) FROM taxonomy").fetchone()[0]
            cache = connection.execute("SELECT COUNT(*) FROM genbank_cache").fetchone()[0]
        return {"taxonomy": int(taxonomy), "genbank_cache": int(cache)}

    def marker_counts(self, query: TaxonQuery | None = None) -> dict[str, int]:
        where = ""
        params: list[Any] = []
        if query is not None:
            query_where, params = _query_filters(query, "t")
            where = f"WHERE {' AND '.join(query_where)}"
        sql = f"""
        SELECT
            SUM(c.has_its) AS its,
            SUM(c.has_its2) AS its2,
            SUM(c.has_matk) AS matk,
            SUM(c.has_rbcl) AS rbcl
        FROM genbank_cache c
        JOIN taxonomy t ON t.taxon_id = c.taxon_id
        {where}
        """
        with self.connect() as connection:
            row = connection.execute(sql, params).fetchone()
        return {key: int(row[key] or 0) for key in row.keys()}

    def taxon_summaries(
        self,
        rank: str,
        query: TaxonQuery | None = None,
    ) -> list[dict[str, Any]]:
        if rank not in {"family", "genus", "species"}:
            raise ValueError(f"unsupported taxon rank: {rank}")
        name_column = "scientific_name" if rank == "species" else rank
        where = [f"t.{name_column} IS NOT NULL", f"trim(t.{name_column}) != ''"]
        params: list[Any] = []
        if query is not None:
            query_where, query_params = _query_filters(query, "t")
            where.extend(query_where)
            params.extend(query_params)
        sql = f"""
        SELECT
            MIN(t.{name_column}) AS name,
            COUNT(*) AS records,
            SUM(c.has_its) AS its,
            SUM(c.has_its2) AS its2,
            SUM(c.has_matk) AS matk,
            SUM(c.has_rbcl) AS rbcl
        FROM genbank_cache c
        JOIN taxonomy t ON t.taxon_id = c.taxon_id
        WHERE {' AND '.join(where)}
        GROUP BY lower(t.{name_column})
        ORDER BY lower(name)
        """
        with self.connect() as connection:
            rows = connection.execute(sql, params).fetchall()
        return [
            {
                "name": row["name"],
                "records": int(row["records"]),
                "markers": {
                    "its": int(row["its"] or 0),
                    "its2": int(row["its2"] or 0),
                    "matk": int(row["matk"] or 0),
                    "rbcl": int(row["rbcl"] or 0),
                },
            }
            for row in rows
        ]


def _query_filters(query: TaxonQuery, alias: str) -> tuple[list[str], list[Any]]:
    filters = [(query.rank, query.name)]
    filters.extend(
        (constraint.rank, constraint.name)
        for constraint in getattr(query, "constraints", ())
    )

    where: list[str] = []
    params: list[Any] = []
    seen: set[tuple[str, str]] = set()
    for rank, name in filters:
        column = TAXONOMY_FILTER_COLUMNS.get(rank)
        if column is None:
            raise ValueError(f"unsupported taxon rank: {rank}")
        key = (rank, name.casefold())
        if key in seen:
            continue
        seen.add(key)
        where.append(f"lower({alias}.{column}) = lower(?)")
        params.append(name)
    return where, params


def _cache_from_row(row: sqlite3.Row) -> GenBankCacheRecord:
    updated_at = row["updated_at"]
    return GenBankCacheRecord(
        accession_root=row["accession_root"],
        version=int(row["version"]),
        accession_version=row["accession_version"],
        taxon_id=int(row["taxon_id"]),
        has_its=bool(row["has_its"]),
        has_matk=bool(row["has_matk"]),
        has_rbcl=bool(row["has_rbcl"]),
        has_its2=bool(row["has_its2"]),
        updated_at=datetime.fromisoformat(updated_at) if updated_at else None,
    )


def _taxonomy_from_joined_row(row: sqlite3.Row) -> TaxonomyRecord:
    return TaxonomyRecord(
        taxon_id=int(row["t_taxon_id"]),
        scientific_name=row["scientific_name"],
        kingdom=row["kingdom"],
        phylum=row["phylum"],
        class_name=row["class"],
        order=row["order"],
        family=row["family"],
        genus=row["genus"],
        species=row["species"],
        infraspecific_rank=row["infraspecific_rank"],
        is_hybrid=bool(row["is_hybrid"]),
        is_uncertain=bool(row["is_uncertain"]),
    )


def _chunks(values: list[T], size: int) -> Iterable[list[T]]:
    for index in range(0, len(values), size):
        yield values[index : index + size]
