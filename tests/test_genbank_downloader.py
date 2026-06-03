from __future__ import annotations

import inspect
from pathlib import Path
from typing import Any

import pytest
from requests.adapters import HTTPAdapter

import barcode_kit.genbank as genbank_module
from barcode_kit.exceptions import GenBankError, GenBankFetchError, GenBankSearchError
from barcode_kit.genbank import NCBIGenBankClient


class FakeSearchResponse:
    def __init__(self, text: str, *, status_code: int = 200):
        self.text = text
        self.status_code = status_code


class FakeSearchSession:
    def __init__(self, responses: list[FakeSearchResponse]):
        self.responses = responses
        self.requests: list[tuple[str, str]] = []

    def request(self, method: str, url: str, timeout: float, **kwargs):
        self.requests.append((method, url))
        return self.responses.pop(0)

    def close(self) -> None:
        pass


class FakeStreamResponse:
    def __init__(self, text: str, *, status_code: int = 200, chunk_size: int = 17):
        self.text = text
        self.status_code = status_code
        self.chunk_size = chunk_size

    def iter_content(self, chunk_size: int, decode_unicode: bool = False):
        data = self.text.encode("utf-8")
        for index in range(0, len(data), self.chunk_size):
            yield data[index : index + self.chunk_size]


class FakeDownloadSession:
    def __init__(self, responses: list[FakeStreamResponse]):
        self.responses = responses
        self.posted_ids: list[str] = []

    def post(self, url: str, data: dict[str, Any], stream: bool, timeout: float):
        self.posted_ids.append(str(data["id"]))
        return self.responses.pop(0)

    def close(self) -> None:
        pass


def build_client_with_fake_sessions(
    monkeypatch,
    download_session: FakeDownloadSession,
    **kwargs,
) -> NCBIGenBankClient:
    search_session = genbank_module._build_ncbi_session()
    sessions = [search_session, download_session]
    monkeypatch.setattr(genbank_module, "_build_ncbi_session", lambda: sessions.pop(0))
    return NCBIGenBankClient(
        email="test@example.com",
        download_workers=1,
        **kwargs,
    )


def test_genbank_module_does_not_export_callback_aliases():
    assert "DownloadStatus" in genbank_module.__all__
    assert "DownloadProgressCallback" not in genbank_module.__all__
    assert "DownloadedRecordCallback" not in genbank_module.__all__


def test_search_does_not_retry_non_rate_limited_client_errors(monkeypatch):
    session = FakeSearchSession([FakeSearchResponse("bad request", status_code=400)])
    monkeypatch.setattr(genbank_module, "_build_ncbi_session", lambda: session)
    client = NCBIGenBankClient(email="test@example.com", retry_attempts=3)

    try:
        with pytest.raises(GenBankSearchError, match="HTTP 400") as exc_info:
            client.search_accessions("bad term")

        assert isinstance(exc_info.value, GenBankError)
        assert len(session.requests) == 1
    finally:
        client.close()


def test_search_retries_use_configured_retry_attempts(monkeypatch):
    session = FakeSearchSession(
        [
            FakeSearchResponse("temporary", status_code=500),
            FakeSearchResponse("temporary", status_code=500),
        ]
    )
    monkeypatch.setattr(genbank_module, "_build_ncbi_session", lambda: session)
    client = NCBIGenBankClient(email="test@example.com", retry_attempts=2)

    try:
        with pytest.raises(GenBankSearchError, match="HTTP 500"):
            client.search_accessions("Iris[Organism] AND rbcl")

        assert len(session.requests) == 2
    finally:
        client.close()


def test_fetch_response_framing_errors_use_fetch_exception():
    response = FakeStreamResponse("not a GenBank record\n")

    with pytest.raises(GenBankFetchError, match="unexpected response body") as exc_info:
        list(genbank_module._iter_genbank_record_texts(response))

    assert isinstance(exc_info.value, GenBankError)


def test_fetch_records_streams_complete_records_to_cache(monkeypatch, tmp_path: Path, genbank_text):
    first = genbank_text(accession="PP476489", version=4)
    second = genbank_text(accession="QQ123456", version=1)
    session = FakeDownloadSession([FakeStreamResponse(first + second)])
    client = build_client_with_fake_sessions(
        monkeypatch,
        session,
        batch_size=2,
        retry_attempts=1,
    )

    report = client.fetch_records(["PP476489.4", "QQ123456.1"], tmp_path)

    assert set(report.succeeded) == {"PP476489.4", "QQ123456.1"}
    assert report.failed == {}
    for accession, record in report.succeeded.items():
        assert record.accession == accession
        assert record.status == "succeeded"
        assert (tmp_path / f"{accession}.gb").read_text(encoding="utf-8").startswith("LOCUS")


def test_fetch_records_reports_streamed_bytes_and_completed_records(
    monkeypatch,
    tmp_path: Path,
    genbank_text,
):
    record = genbank_text(accession="PP476489", version=4)
    byte_events: list[int] = []
    completed_accessions: list[str] = []
    session = FakeDownloadSession(
        [FakeStreamResponse(record, chunk_size=max(1, len(record.encode("utf-8")) // 3))]
    )
    client = build_client_with_fake_sessions(
        monkeypatch,
        session,
        batch_size=1,
        retry_attempts=1,
    )

    report = client.fetch_records(
        ["PP476489.4"],
        tmp_path,
        progress_callback=byte_events.append,
        record_callback=completed_accessions.append,
    )

    assert set(report.succeeded) == {"PP476489.4"}
    assert sum(byte_events) == len(record.encode("utf-8"))
    assert len(byte_events) > 1
    assert completed_accessions == ["PP476489.4"]


def test_fetch_records_retries_only_missing_accessions(monkeypatch, tmp_path: Path, genbank_text):
    first = genbank_text(accession="PP476489", version=4)
    second = genbank_text(accession="QQ123456", version=1)
    session = FakeDownloadSession(
        [
            FakeStreamResponse(first),
            FakeStreamResponse(second),
        ]
    )
    client = build_client_with_fake_sessions(
        monkeypatch,
        session,
        batch_size=2,
        retry_attempts=2,
    )

    report = client.fetch_records(["PP476489.4", "QQ123456.1"], tmp_path)

    assert set(report.succeeded) == {"PP476489.4", "QQ123456.1"}
    assert report.failed == {}
    assert session.posted_ids == ["PP476489.4,QQ123456.1", "QQ123456.1"]


@pytest.mark.parametrize(
    "response_text",
    [
        "",
        "LOCUS       PP476489\n",
    ],
)
def test_fetch_records_retries_transient_response_framing_errors(
    monkeypatch,
    tmp_path: Path,
    genbank_text,
    response_text: str,
):
    record = genbank_text(accession="PP476489", version=4)
    session = FakeDownloadSession(
        [
            FakeStreamResponse(response_text),
            FakeStreamResponse(record),
        ]
    )
    client = build_client_with_fake_sessions(
        monkeypatch,
        session,
        batch_size=1,
        retry_attempts=2,
    )

    report = client.fetch_records(["PP476489.4"], tmp_path)

    assert set(report.succeeded) == {"PP476489.4"}
    assert report.failed == {}
    assert session.posted_ids == ["PP476489.4", "PP476489.4"]


def test_fetch_records_reports_unified_items_with_attempt_counts(
    monkeypatch,
    tmp_path: Path,
    genbank_text,
):
    first = genbank_text(accession="PP476489", version=4)
    session = FakeDownloadSession([FakeStreamResponse(first)])
    client = build_client_with_fake_sessions(
        monkeypatch,
        session,
        batch_size=2,
        retry_attempts=1,
    )

    report = client.fetch_records(["PP476489.4", "QQ123456.1"], tmp_path)

    items = getattr(report, "items", None)
    assert items is not None
    assert set(items) == {"PP476489.4", "QQ123456.1"}
    assert items["PP476489.4"].status == "succeeded"
    assert items["PP476489.4"].attempts == 1
    assert not hasattr(items["PP476489.4"], "path")
    assert not hasattr(items["PP476489.4"], "size")
    assert items["QQ123456.1"].status == "failed"
    assert items["QQ123456.1"].error_type == "missing_from_response"
    assert items["QQ123456.1"].retryable is True
    assert items["QQ123456.1"].attempts == 1


def test_fetch_records_reuses_valid_cached_file(monkeypatch, tmp_path: Path, genbank_text):
    cached_path = tmp_path / "PP476489.4.gb"
    cached_path.write_text(genbank_text(accession="PP476489", version=4), encoding="utf-8")
    session = FakeDownloadSession([])
    client = build_client_with_fake_sessions(
        monkeypatch,
        session,
        batch_size=1,
        retry_attempts=1,
    )

    report = client.fetch_records(["PP476489.4"], tmp_path)

    assert set(report.succeeded) == {"PP476489.4"}
    assert report.succeeded["PP476489.4"].status == "succeeded"
    assert cached_path.exists()
    assert session.posted_ids == []


def test_fetch_records_uses_deque_for_pending_work():
    source = inspect.getsource(NCBIGenBankClient.fetch_records)

    assert "deque[DownloadItem]" in source
    assert "popleft()" in source


def test_client_exposes_only_fetch_api_for_record_downloads():
    removed_name = "download" + "_records"
    assert not hasattr(NCBIGenBankClient, removed_name)


def test_client_does_not_accept_external_search_session():
    signature = inspect.signature(NCBIGenBankClient)

    assert "session" not in signature.parameters
    assert "worker_session_factory" not in signature.parameters


def test_client_configures_internal_search_session_with_single_connection_pool():
    client = NCBIGenBankClient(email="test@example.com")

    adapter = client._session.get_adapter("https://")

    assert isinstance(adapter, HTTPAdapter)
    assert adapter._pool_connections == 1
    assert adapter._pool_maxsize == 1
    assert adapter._pool_block is True
    assert adapter.max_retries.total == 0
