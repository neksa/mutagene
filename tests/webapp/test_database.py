"""Tests for webapp database models and manager."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from mutagene.webapp.database.manager import DatabaseManager
from mutagene.webapp.database.models import Analysis, File, Result, _SafeEncoder, init_db


@pytest.fixture
def db(tmp_path):
    """Create a temporary database."""
    return DatabaseManager(db_path=tmp_path / "test.db")


class TestSafeEncoder:
    def test_numpy_int(self):
        assert json.dumps({"x": np.int64(42)}, cls=_SafeEncoder) == '{"x": 42}'

    def test_numpy_float(self):
        result = json.loads(json.dumps({"x": np.float64(3.14)}, cls=_SafeEncoder))
        assert abs(result["x"] - 3.14) < 1e-10

    def test_numpy_array(self):
        result = json.loads(json.dumps({"x": np.array([1, 2, 3])}, cls=_SafeEncoder))
        assert result["x"] == [1, 2, 3]

    def test_plain_types_pass_through(self):
        data = {"a": 1, "b": "hello", "c": [1, 2], "d": None}
        assert json.loads(json.dumps(data, cls=_SafeEncoder)) == data


class TestInitDb:
    def test_creates_tables(self, tmp_path):
        db_path = tmp_path / "test.db"
        init_db(str(db_path))
        assert db_path.exists()

        import sqlite3

        conn = sqlite3.connect(str(db_path))
        tables = [
            r[0]
            for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
        ]
        conn.close()
        assert "analyses" in tables
        assert "results" in tables
        assert "files" in tables

    def test_idempotent(self, tmp_path):
        db_path = tmp_path / "test.db"
        init_db(str(db_path))
        init_db(str(db_path))  # should not raise


class TestAnalysisCRUD:
    def test_create_and_get(self, db):
        aid = db.create_analysis("test.maf", "hg19")
        analysis = db.get_analysis(aid)
        assert analysis["name"] == "test.maf"
        assert analysis["genome"] == "hg19"
        assert analysis["signatures"] == "COSMICv3"
        assert analysis["status"] == "pending"

    def test_create_with_config(self, db):
        config = {"classify": True, "cluster": False}
        aid = db.create_analysis("test.maf", "hg38", config=config)
        analysis = db.get_analysis(aid)
        assert analysis["config"] == config

    def test_get_nonexistent(self, db):
        assert db.get_analysis(999) is None

    def test_update_status(self, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.update_analysis_status(aid, "running")
        assert db.get_analysis(aid)["status"] == "running"

    def test_update_status_with_error(self, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.update_analysis_status(aid, "error", "something broke")
        analysis = db.get_analysis(aid)
        assert analysis["status"] == "error"
        assert analysis["error_message"] == "something broke"

    def test_update_counts(self, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.update_analysis_counts(aid, samples=3, mutations=150)
        analysis = db.get_analysis(aid)
        assert analysis["samples"] == 3
        assert analysis["mutations"] == 150

    def test_list_analyses(self, db):
        db.create_analysis("a.maf", "hg19")
        db.create_analysis("b.maf", "hg38")
        db.create_analysis("c.maf", "hg19")
        analyses = db.list_analyses()
        assert len(analyses) == 3
        names = {a["name"] for a in analyses}
        assert names == {"a.maf", "b.maf", "c.maf"}

    def test_list_analyses_limit(self, db):
        for i in range(5):
            db.create_analysis(f"{i}.maf", "hg19")
        assert len(db.list_analyses(limit=2)) == 2

    def test_delete(self, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.delete_analysis(aid)
        assert db.get_analysis(aid) is None


class TestResultCRUD:
    def test_store_and_retrieve(self, db):
        aid = db.create_analysis("test.maf", "hg19")
        data = {"signatures": {"SBS1": 0.3, "SBS5": 0.7}}
        db.store_result(aid, "signatures", data)
        results = db.get_results_by_type(aid, "signatures")
        assert len(results) == 1
        assert results[0]["data"] == data

    def test_store_numpy_data(self, db):
        aid = db.create_analysis("test.maf", "hg19")
        data = {"values": np.array([0.1, 0.2, 0.3]), "count": np.int64(42)}
        db.store_result(aid, "profile", data)
        results = db.get_results_by_type(aid, "profile")
        assert results[0]["data"]["values"] == [0.1, 0.2, 0.3]
        assert results[0]["data"]["count"] == 42

    def test_get_all_results(self, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.store_result(aid, "profile", {"p": 1})
        db.store_result(aid, "signatures", {"s": 2})
        results = db.get_all_results(aid)
        assert len(results) == 2

    def test_empty_results(self, db):
        aid = db.create_analysis("test.maf", "hg19")
        assert db.get_results_by_type(aid, "profile") == []
        assert db.get_all_results(aid) == []


class TestFileCRUD:
    def test_register_and_get(self, db, tmp_path):
        aid = db.create_analysis("test.maf", "hg19")
        test_file = tmp_path / "input.maf"
        test_file.write_text("data")
        fid = db.register_file(aid, "input_maf", "input.maf", str(test_file))
        record = db.get_file(fid)
        assert record["filename"] == "input.maf"
        assert record["file_type"] == "input_maf"
        assert record["size_bytes"] == 4

    def test_get_by_type(self, db, tmp_path):
        aid = db.create_analysis("test.maf", "hg19")
        f1 = tmp_path / "a.maf"
        f1.write_text("x")
        f2 = tmp_path / "b.tsv"
        f2.write_text("y")
        db.register_file(aid, "input_maf", "a.maf", str(f1))
        db.register_file(aid, "profile_tsv", "b.tsv", str(f2))
        maf_files = db.get_files_by_type(aid, "input_maf")
        assert len(maf_files) == 1
        assert maf_files[0]["filename"] == "a.maf"

    def test_get_all_files(self, db, tmp_path):
        aid = db.create_analysis("test.maf", "hg19")
        for name in ["a.maf", "b.tsv", "c.png"]:
            f = tmp_path / name
            f.write_text("x")
            db.register_file(aid, "output", name, str(f))
        assert len(db.get_all_files(aid)) == 3

    def test_get_nonexistent_file(self, db):
        assert db.get_file(999) is None
