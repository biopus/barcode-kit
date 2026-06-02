from __future__ import annotations

import codecs
import json
import re
import threading
import time
from collections import deque
from collections.abc import Callable, Iterable
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

import requests
from requests.adapters import HTTPAdapter
from requests.exceptions import ChunkedEncodingError, ContentDecodingError, ReadTimeout

from barcode_kit.config import AppConfig, ensure_app_dirs
from barcode_kit.exceptions import (
    ConfigError,
    GenBankError,
    GenBankFetchError,
    GenBankSearchError,
    TaxonomyError,
)
from barcode_kit.models import GenBankCacheRecord, Marker, TaxonQuery
from barcode_kit.parser import (
    AccessionVersion,
    parse_accession_version,
    parse_genbank_file,
    parse_genbank_text,
)
from barcode_kit.storage import Storage


__all__ = [
    "DownloadItem",
    "DownloadProgressCallback",
    "DownloadReport",
    "DownloadStatus",
    "DownloadedRecordCallback",
    "NCBIGenBankClient",
    "RequestLimiter",
    "SyncResult",
    "SyncService",
]


BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
FILENAME_RE = re.compile(r"[^A-Za-z0-9._-]+")
DownloadStatus = Literal["pending", "succeeded", "failed"]
DownloadProgressCallback = Callable[[int], None]
DownloadedRecordCallback = Callable[[str], None]


def _build_ncbi_session() -> requests.Session:
    session = requests.Session()
    adapter = HTTPAdapter(
        pool_connections=1,
        pool_maxsize=1,
        pool_block=True,
        max_retries=0,
    )
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


@dataclass(frozen=True)
class SyncResult:
    query: str
    remote_count: int
    downloaded: int
    reused_local: int
    ingested: int
    skipped: int
    updated: int
    failed: list[dict[str, str]]


@dataclass
class DownloadItem:
    accession: str
    attempts: int = 0
    status: DownloadStatus = "pending"
    error_type: str | None = None
    message: str | None = None
    retryable: bool = False

    @property
    def succeeded(self) -> bool:
        return self.status == "succeeded"

    @property
    def failed(self) -> bool:
        return self.status == "failed"

    def mark_succeeded(self) -> None:
        self.status = "succeeded"
        self.error_type = None
        self.message = None
        self.retryable = False

    def mark_failed(self, error_type: str, message: str, retryable: bool) -> None:
        self.status = "failed"
        self.error_type = error_type
        self.message = message
        self.retryable = retryable


@dataclass(frozen=True)
class DownloadReport:
    items: dict[str, DownloadItem]

    @property
    def succeeded(self) -> dict[str, DownloadItem]:
        return {accession: item for accession, item in self.items.items() if item.succeeded}

    @property
    def failed(self) -> dict[str, DownloadItem]:
        return {accession: item for accession, item in self.items.items() if item.failed}


class RequestLimiter:
    def __init__(self, min_interval: float):
        self.min_interval = min_interval
        self._lock = threading.Lock()
        self._last_started_at = 0.0

    def wait_turn(self) -> None:
        if self.min_interval <= 0:
            return
        with self._lock:
            elapsed = time.monotonic() - self._last_started_at
            if elapsed < self.min_interval:
                time.sleep(self.min_interval - elapsed)
            self._last_started_at = time.monotonic()


class NCBIGenBankClient:
    def __init__(
        self,
        email: str,
        api_key: str | None = None,
        *,
        batch_size: int = 500,
        timeout: float = 30,
        retry_attempts: int = 3,
        base_url: str = BASE_URL,
        download_workers: int = 1,
    ):
        if not email:
            raise ConfigError("collectors.genbank.email is required for NCBI E-utilities sync")
        self.email = email
        self.api_key = api_key
        self.batch_size = batch_size
        self.timeout = timeout
        self.retry_attempts = retry_attempts
        self.base_url = base_url
        self._session = _build_ncbi_session()
        self._min_interval = 0.12 if api_key else 0.38
        self._limiter = RequestLimiter(self._min_interval)
        self.download_workers = max(1, download_workers)

    def close(self) -> None:
        self._session.close()

    def search_accessions(self, term: str) -> list[str]:
        accessions: list[str] = []
        retstart = 0
        while True:
            params = {
                **self._common_params(),
                "term": term,
                "retmode": "json",
                "idtype": "acc",
                "retstart": retstart,
                "retmax": 10000,
            }
            data = self._request_json("GET", "esearch.fcgi", params=params)
            result = data.get("esearchresult", {})
            count = int(result.get("count", 0))
            idlist = [str(value) for value in result.get("idlist", [])]
            accessions.extend(idlist)
            if count == 0 or len(accessions) >= count or not idlist:
                return _dedupe(accessions)
            retstart += 10000

    def fetch_records(
        self,
        accessions: list[str],
        output_dir: Path,
        *,
        progress_callback: DownloadProgressCallback | None = None,
        record_callback: DownloadedRecordCallback | None = None,
    ) -> DownloadReport:
        output_dir.mkdir(parents=True, exist_ok=True)
        items: dict[str, DownloadItem] = {}
        pending: deque[DownloadItem] = deque()
        for accession in _dedupe(accessions):
            cached = _valid_cached_record(output_dir, accession)
            if cached is not None:
                items[accession] = cached
                continue
            item = DownloadItem(accession=accession)
            items[accession] = item
            pending.append(item)

        if not pending:
            return DownloadReport(items=items)

        thread_local = threading.local()
        created_sessions: list[Any] = []
        created_sessions_lock = threading.Lock()

        def get_worker_session() -> requests.Session | Any:
            session = getattr(thread_local, "session", None)
            if session is not None:
                return session
            session = _build_ncbi_session()
            thread_local.session = session
            with created_sessions_lock:
                created_sessions.append(session)
            return session

        def run_batch(batch: list[DownloadItem]) -> DownloadReport:
            return self._download_batch(
                batch,
                output_dir,
                get_worker_session(),
                progress_callback=progress_callback,
                record_callback=record_callback,
            )

        in_flight: dict[Any, list[DownloadItem]] = {}
        try:
            with ThreadPoolExecutor(max_workers=self.download_workers) as executor:
                while pending or in_flight:
                    while pending and len(in_flight) < self.download_workers:
                        batch_items = [
                            pending.popleft()
                            for _ in range(min(self.batch_size, len(pending)))
                        ]
                        for item in batch_items:
                            item.attempts += 1
                        in_flight[executor.submit(run_batch, batch_items)] = batch_items

                    done, _ = wait(in_flight, return_when=FIRST_COMPLETED)
                    for future in done:
                        in_flight.pop(future)
                        result = future.result()
                        items.update(result.items)
                        for item in result.failed.values():
                            if item.retryable and item.attempts < self.retry_attempts:
                                pending.append(item)
        finally:
            for session in created_sessions:
                close = getattr(session, "close", None)
                if callable(close):
                    close()

        return DownloadReport(items=items)

    def _download_batch(
        self,
        batch: list[DownloadItem],
        output_dir: Path,
        session: requests.Session | Any,
        *,
        progress_callback: DownloadProgressCallback | None = None,
        record_callback: DownloadedRecordCallback | None = None,
    ) -> DownloadReport:
        requested = {item.accession for item in batch}
        item_by_accession = {item.accession: item for item in batch}
        data = {
            **self._common_params(),
            "id": ",".join(item.accession for item in batch),
            "rettype": "gbwithparts",
            "retmode": "text",
        }
        self._limiter.wait_turn()
        try:
            response = session.post(
                self.base_url + "efetch.fcgi",
                data=data,
                stream=True,
                timeout=self.timeout,
            )
        except ReadTimeout as error:
            return _failed_batch(batch, "read_timeout", str(error), True)
        except requests.ConnectionError as error:
            return _failed_batch(batch, "connection_error", str(error), True)
        except requests.RequestException as error:
            return _failed_batch(batch, type(error).__name__, str(error), True)

        status_failure = _http_failure(response)
        if status_failure is not None:
            error_type, message, retryable = status_failure
            return _failed_batch(
                batch,
                error_type,
                message,
                retryable,
            )

        try:
            for record_text in _iter_genbank_record_texts(
                response,
                progress_callback=progress_callback,
            ):
                parsed = parse_genbank_text(record_text)
                accession = parsed.accession.value
                if accession not in requested:
                    return _failed_batch(
                        batch,
                        "accession_mismatch",
                        f"record accession {accession} not present in current batch",
                        False,
                    )
                item = item_by_accession[accession]
                if item.succeeded:
                    continue
                _write_downloaded_record(output_dir, accession, record_text)
                item.mark_succeeded()
                if record_callback is not None:
                    record_callback(accession)
        except GenBankFetchError as error:
            return _failed_batch(
                batch,
                error.error_type,
                str(error),
                error.retryable,
            )
        except GenBankError as error:
            return _failed_batch(
                batch,
                "parse_error",
                str(error),
                False,
            )
        except UnicodeDecodeError as error:
            return _failed_batch(
                batch,
                "decode_error",
                str(error),
                True,
            )
        except (ChunkedEncodingError, ContentDecodingError, requests.ConnectionError) as error:
            return _failed_batch(
                batch,
                "stream_incomplete",
                str(error),
                True,
            )
        except OSError as error:
            return _failed_batch(
                batch,
                "write_error",
                str(error),
                False,
            )

        for item in batch:
            if not item.succeeded:
                item.mark_failed(
                    "missing_from_response",
                    "EFetch response did not include requested accession",
                    True,
                )
        return DownloadReport(items={item.accession: item for item in batch})

    def _request_json(self, method: str, endpoint: str, **kwargs) -> dict:
        text = self._request_text(method, endpoint, **kwargs)
        try:
            return json.loads(text)
        except json.JSONDecodeError as error:
            raise GenBankSearchError(f"NCBI returned invalid JSON from {endpoint}") from error

    def _request_text(self, method: str, endpoint: str, **kwargs) -> str:
        last_error: Exception | None = None
        for attempt in range(1, max(1, self.retry_attempts) + 1):
            self._limiter.wait_turn()
            try:
                response = self._session.request(
                    method,
                    self.base_url + endpoint,
                    timeout=self.timeout,
                    **kwargs,
                )
            except requests.RequestException as error:
                last_error = error
                if attempt >= max(1, self.retry_attempts):
                    break
                continue

            status_failure = _http_failure(response)
            if status_failure is None:
                return response.text

            _, message, retryable = status_failure
            error = GenBankSearchError(message)
            if not retryable:
                raise error
            last_error = error
            if attempt >= max(1, self.retry_attempts):
                break
        raise GenBankSearchError(str(last_error) if last_error else "NCBI request failed") from last_error

    def _common_params(self) -> dict[str, str]:
        params = {"db": "nuccore", "tool": "barcode-kit", "email": self.email}
        if self.api_key:
            params["api_key"] = self.api_key
        return params


class SyncService:
    def __init__(
        self,
        config: AppConfig,
        storage: Storage,
        taxonomy_resolver: Any,
        client: Any | None = None,
    ):
        self.config = config
        self.storage = storage
        self.taxonomy_resolver = taxonomy_resolver
        self.client = client or NCBIGenBankClient(
            email=config.collectors.genbank_email,
            api_key=config.collectors.genbank_api_key,
            batch_size=config.collectors.batch_size,
            timeout=config.collectors.timeout,
            retry_attempts=config.collectors.retry_attempts,
            download_workers=config.collectors.download_workers,
        )

    def sync(self, query: TaxonQuery, marker: Marker, progress: Any | None = None) -> SyncResult:
        ensure_app_dirs(self.config)
        self.storage.initialize()
        term = f"{query.ncbi_term()} AND {marker.ncbi_search_term}"
        remote_accessions = [parse_accession_version(value) for value in self.client.search_accessions(term)]
        local_versions = self.storage.get_cached_versions_by_roots(
            accession.root for accession in remote_accessions
        )
        to_download: list[AccessionVersion] = []
        skipped = 0
        updated = 0
        for accession in remote_accessions:
            cached_version = local_versions.get(accession.root)
            if cached_version is None:
                to_download.append(accession)
            elif accession.version > cached_version:
                to_download.append(accession)
                updated += 1
            else:
                skipped += 1

        failed: list[dict[str, str]] = []
        reused_local = 0
        downloaded = 0
        ingested = 0
        if to_download:
            record_accessions: list[str] = []
            requested_remote: list[str] = []
            for accession in to_download:
                cached = _valid_cached_record(self.config.genbank_cache_dir, accession.value)
                if cached is None:
                    requested_remote.append(accession.value)
                    continue
                record_accessions.append(accession.value)
                reused_local += 1

            if requested_remote:
                if progress is not None:
                    progress.start_download(total_records=len(requested_remote))
                    report = self.client.fetch_records(
                        requested_remote,
                        self.config.genbank_cache_dir,
                        progress_callback=progress.record_downloaded_bytes,
                        record_callback=progress.record_downloaded_record,
                    )
                    progress.finish_download(
                        succeeded_records=len(report.succeeded),
                        failed_records=len(report.failed),
                    )
                else:
                    report = self.client.fetch_records(
                        requested_remote,
                        self.config.genbank_cache_dir,
                    )
                record_accessions.extend(report.succeeded)
                downloaded += len(report.succeeded)
                failed.extend(
                    {
                        "accession": failure.accession,
                        "stage": "download",
                        "error": f"{failure.error_type}: {failure.message}",
                    }
                    for failure in report.failed.values()
                )

            for accession_value in record_accessions:
                try:
                    self._ingest_record_file(
                        accession_value,
                        _cached_record_path(self.config.genbank_cache_dir, accession_value),
                    )
                    ingested += 1
                except (GenBankError, TaxonomyError) as error:
                    failed.append(
                        {"accession": accession_value, "stage": "ingest", "error": str(error)}
                    )

        result = SyncResult(
            query=term,
            remote_count=len(remote_accessions),
            downloaded=downloaded,
            reused_local=reused_local,
            ingested=ingested,
            skipped=skipped,
            updated=updated,
            failed=failed,
        )
        self._write_sync_log(result)
        return result

    def close(self) -> None:
        close = getattr(self.client, "close", None)
        if callable(close):
            close()

    def _ingest_record_file(self, requested_accession: str, path: Path) -> None:
        parsed = parse_genbank_file(path)
        if parsed.accession.value != requested_accession:
            raise GenBankError(
                f"requested {requested_accession} but parsed {parsed.accession.value}"
            )
        taxonomy = self.taxonomy_resolver.standardize(parsed.organism, parsed.taxon_id)
        cache_record = GenBankCacheRecord(
            accession_root=parsed.accession.root,
            version=parsed.accession.version,
            accession_version=parsed.accession.value,
            taxon_id=taxonomy.taxon_id,
            has_its=parsed.marker_flags[Marker.ITS],
            has_matk=parsed.marker_flags[Marker.MATK],
            has_rbcl=parsed.marker_flags[Marker.RBCL],
            has_its2=parsed.marker_flags[Marker.ITS2],
            updated_at=datetime.now(timezone.utc),
        )
        with self.storage.connect() as connection:
            self.storage.upsert_taxonomy(taxonomy, connection)
            self.storage.upsert_genbank_cache(cache_record, connection)

    def _write_sync_log(self, result: SyncResult) -> None:
        self.config.logs_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        path = self.config.logs_dir / f"sync-{stamp}.json"
        path.write_text(
            json.dumps(asdict(result), ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )


def _iter_genbank_record_texts(
    response: Any,
    chunk_size: int = 64 * 1024,
    *,
    progress_callback: DownloadProgressCallback | None = None,
) -> Iterable[str]:
    current: list[str] = []
    in_record = False
    saw_any_line = False
    decoder = codecs.getincrementaldecoder("utf-8")()
    pending_text = ""

    def process_line(line: str) -> str | None:
        nonlocal current, in_record, saw_any_line
        saw_any_line = True
        if line.startswith("LOCUS"):
            if in_record:
                raise GenBankFetchError("record started before previous terminator")
            current = [line]
            in_record = True
            return None
        if in_record:
            current.append(line)
            if line.strip() == "//":
                record_text = "\n".join(current) + "\n"
                current = []
                in_record = False
                return record_text
            return None
        if line.strip():
            raise GenBankFetchError("unexpected response body")
        return None

    for chunk in response.iter_content(chunk_size=chunk_size, decode_unicode=False):
        if not chunk:
            continue
        if progress_callback is not None:
            progress_callback(len(chunk))
        pending_text += decoder.decode(chunk)
        lines = pending_text.splitlines(keepends=True)
        pending_text = ""
        if lines and not lines[-1].endswith(("\n", "\r")):
            pending_text = lines.pop()
        for raw_line in lines:
            record_text = process_line(raw_line.rstrip("\r\n"))
            if record_text is not None:
                yield record_text

    pending_text += decoder.decode(b"", final=True)
    if pending_text:
        record_text = process_line(pending_text.rstrip("\r\n"))
        if record_text is not None:
            yield record_text
    if not saw_any_line:
        raise GenBankFetchError("empty response body", error_type="empty_response", retryable=True)
    if in_record:
        raise GenBankFetchError(
            "stream ended before record terminator",
            error_type="stream_incomplete",
            retryable=True,
        )


def _write_downloaded_record(output_dir: Path, accession: str, record_text: str) -> Path:
    parsed = parse_genbank_text(record_text)
    if parsed.accession.value != accession:
        raise GenBankFetchError(
            f"requested {accession} but parsed {parsed.accession.value}",
            error_type="accession_mismatch",
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = _cached_record_path(output_dir, accession)
    tmp_path = output_dir / f"{_safe_filename(accession)}.gb.tmp-{threading.get_ident()}"
    try:
        tmp_path.write_text(record_text, encoding="utf-8")
        parsed_tmp = parse_genbank_file(tmp_path)
        if parsed_tmp.accession.value != accession:
            raise GenBankFetchError(
                f"requested {accession} but parsed {parsed_tmp.accession.value}",
                error_type="accession_mismatch",
            )
        tmp_path.replace(output_path)
    finally:
        tmp_path.unlink(missing_ok=True)
    return output_path


def _valid_cached_record(output_dir: Path, accession: str) -> DownloadItem | None:
    path = _cached_record_path(output_dir, accession)
    if not path.exists():
        return None
    try:
        parsed = parse_genbank_file(path)
    except GenBankError:
        path.unlink(missing_ok=True)
        return None
    if parsed.accession.value != accession:
        path.unlink(missing_ok=True)
        return None
    return DownloadItem(accession=accession, status="succeeded")


def _cached_record_path(output_dir: Path, accession: str) -> Path:
    return output_dir / f"{_safe_filename(accession)}.gb"


def _http_failure(response: Any) -> tuple[str, str, bool] | None:
    status = int(getattr(response, "status_code", 200))
    if status == 429:
        return ("http_429", "HTTP 429", True)
    if 500 <= status <= 599:
        return ("http_5xx", f"HTTP {status}", True)
    if status >= 400:
        body = str(getattr(response, "text", ""))[:200]
        return (f"http_{status}", f"HTTP {status}: {body}", False)
    return None


def _failed_batch(
    batch: list[DownloadItem],
    error_type: str,
    message: str,
    retryable: bool,
) -> DownloadReport:
    for item in batch:
        if not item.succeeded:
            item.mark_failed(error_type, message, retryable)
    return DownloadReport(items={item.accession: item for item in batch})


def _dedupe(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


def _safe_filename(value: str) -> str:
    return FILENAME_RE.sub("_", value.strip()).strip("._-") or "record"
