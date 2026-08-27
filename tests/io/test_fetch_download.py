"""Resuming a download.

The Range header asked through `file_size`, one byte past the end of an
inclusive range, and the streaming response was never closed. The total size was
also overwritten with the partial response's Content-Length, so a resumed
download reported progress against the size of the remainder.
"""

import io

import pytest

from mutagene.io import fetch


class FakeRaw:
    def __init__(self, payload):
        self.stream = io.BytesIO(payload)

    def read(self, n):
        return self.stream.read(n)


class FakeResponse:
    def __init__(self, payload, status_code=206, headers=None):
        self.status_code = status_code
        self.headers = headers or {"Content-Length": str(len(payload))}
        self.raw = FakeRaw(payload)
        self.closed = False

    def close(self):
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


@pytest.fixture
def captured(monkeypatch):
    """Record the request that download_from_url makes."""
    state = {}

    def fake_head(url):
        return FakeResponse(b"", status_code=200, headers={"Content-Length": "1000"})

    def fake_get(url, headers=None, stream=False):
        state["headers"] = headers
        state["response"] = FakeResponse(b"x" * 400)
        return state["response"]

    monkeypatch.setattr(fetch.requests, "head", fake_head)
    monkeypatch.setattr(fetch.requests, "get", fake_get)
    return state


def test_range_stops_at_the_last_byte(tmp_path, captured):
    dst = tmp_path / "genome.2bit"
    dst.write_bytes(b"y" * 600)  # 600 already downloaded of 1000

    fetch.download_from_url("http://example/genome.2bit", str(dst))

    assert captured["headers"]["Range"] == "bytes=600-999", "requested past the end of the file"


def test_a_fresh_download_starts_at_zero(tmp_path, captured):
    dst = tmp_path / "genome.2bit"

    fetch.download_from_url("http://example/genome.2bit", str(dst))

    assert captured["headers"]["Range"] == "bytes=0-999"


def test_the_response_is_closed(tmp_path, captured):
    dst = tmp_path / "genome.2bit"

    fetch.download_from_url("http://example/genome.2bit", str(dst))

    assert captured["response"].closed, "the connection was never returned to the pool"


def test_the_downloaded_bytes_are_appended(tmp_path, captured):
    dst = tmp_path / "genome.2bit"
    dst.write_bytes(b"y" * 600)

    fetch.download_from_url("http://example/genome.2bit", str(dst))

    assert dst.read_bytes() == b"y" * 600 + b"x" * 400


def test_an_already_complete_file_is_not_refetched(tmp_path, captured):
    dst = tmp_path / "genome.2bit"
    dst.write_bytes(b"y" * 1000)

    fetch.download_from_url("http://example/genome.2bit", str(dst))

    assert "headers" not in captured, "re-requested a file that was already complete"
