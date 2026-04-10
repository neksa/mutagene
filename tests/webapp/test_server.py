"""Tests for webapp Flask server routes."""

import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture
def app(tmp_path):
    """Create test Flask app with temp database and upload folder."""
    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir()
    results_dir = tmp_path / "results"
    results_dir.mkdir()

    with patch("mutagene.webapp.server.GenomeManager") as mock_gm:
        mock_gm_instance = MagicMock()
        mock_gm_instance.get_available_genomes.return_value = ["hg19", "hg38"]
        mock_gm.return_value = mock_gm_instance

        with patch("mutagene.webapp.server.DatabaseManager") as mock_db_cls:
            from mutagene.webapp.database.manager import DatabaseManager

            real_db = DatabaseManager(db_path=tmp_path / "test.db")
            mock_db_cls.return_value = real_db

            from mutagene.webapp.server import create_app

            app, socketio = create_app(
                config={
                    "TESTING": True,
                    "UPLOAD_FOLDER": upload_dir,
                }
            )
            app.config["_DB"] = real_db
            app.config["_UPLOAD_DIR"] = upload_dir
            app.config["_RESULTS_DIR"] = results_dir
            yield app


@pytest.fixture
def client(app):
    return app.test_client()


@pytest.fixture
def db(app):
    return app.config["_DB"]


class TestHomePage:
    def test_index(self, client):
        response = client.get("/")
        assert response.status_code == 200

    def test_history(self, client):
        response = client.get("/history")
        assert response.status_code == 200


class TestUploadAPI:
    def test_upload_maf(self, client, app):
        from io import BytesIO

        data = {
            "file": (BytesIO(b"Hugo_Symbol\tChromosome\n"), "test.maf"),
            "genome": "hg19",
            "signatures": "COSMICv3",
        }
        response = client.post("/api/upload", data=data, content_type="multipart/form-data")
        assert response.status_code == 200
        result = response.get_json()
        assert "analysis_id" in result

    def test_upload_no_file(self, client):
        response = client.post("/api/upload", data={"genome": "hg19"})
        assert response.status_code == 400


class TestAnalysisAPI:
    def test_get_analysis(self, client, db):
        aid = db.create_analysis("test.maf", "hg19")
        response = client.get(f"/api/analysis/{aid}")
        assert response.status_code == 200
        data = response.get_json()
        assert data["name"] == "test.maf"
        assert data["genome"] == "hg19"

    def test_get_analysis_not_found(self, client):
        response = client.get("/api/analysis/999")
        assert response.status_code == 404

    def test_get_analysis_strips_traceback(self, client, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.update_analysis_status(aid, "error", "ValueError: bad input\nTraceback line 1\nline 2")
        response = client.get(f"/api/analysis/{aid}")
        data = response.get_json()
        assert "Traceback" not in data["error_message"]
        assert "bad input" in data["error_message"]


class TestDeleteAPI:
    def test_delete_analysis(self, client, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.update_analysis_status(aid, "complete")
        response = client.delete(f"/api/delete/{aid}")
        assert response.status_code == 200
        assert db.get_analysis(aid) is None

    def test_delete_running_analysis_blocked(self, client, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.update_analysis_status(aid, "running")
        response = client.delete(f"/api/delete/{aid}")
        assert response.status_code == 409


class TestDownloadAPI:
    def test_download_file(self, client, db, app):
        aid = db.create_analysis("test.maf", "hg19")
        upload_dir = app.config["_UPLOAD_DIR"]
        test_file = upload_dir / "output.tsv"
        test_file.write_text("col1\tcol2\n")
        fid = db.register_file(aid, "profile_tsv", "output.tsv", str(test_file))
        response = client.get(f"/api/download/{fid}")
        assert response.status_code == 200

    def test_download_nonexistent(self, client):
        response = client.get("/api/download/999")
        assert response.status_code == 404

    def test_download_outside_allowed_dirs(self, client, db, tmp_path):
        aid = db.create_analysis("test.maf", "hg19")
        # Register a file outside allowed directories
        evil_path = tmp_path / "evil" / "etc_passwd"
        evil_path.parent.mkdir(parents=True)
        evil_path.write_text("root:x:0:0")
        fid = db.register_file(aid, "output", "passwd", str(evil_path))
        response = client.get(f"/api/download/{fid}")
        assert response.status_code == 403


class TestResultsPage:
    def test_results_not_found(self, client):
        response = client.get("/results/999")
        assert response.status_code == 404

    def test_results_analyzing(self, client, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.update_analysis_status(aid, "running")
        response = client.get(f"/results/{aid}")
        assert response.status_code == 200
        assert b"progress" in response.data.lower()
