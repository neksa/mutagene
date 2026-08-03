"""Tests for webapp Flask server routes."""

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


class TestUploadValidation:
    def _upload(self, client, filename="test.maf", **form):
        from io import BytesIO

        data = {"file": (BytesIO(b"Hugo_Symbol\tChromosome\n"), filename)}
        data.update(form)
        return client.post("/api/upload", data=data, content_type="multipart/form-data")

    def test_rejects_unknown_genome(self, client):
        # An unvalidated value here reaches GenomeManager.get_genome_path()
        # and is interpolated into a filesystem path.
        response = self._upload(client, genome="../../../../etc/passwd")
        assert response.status_code == 400
        assert "genome" in response.get_json()["error"].lower()

    def test_rejects_unknown_signature_set(self, client):
        response = self._upload(client, genome="hg19", signatures="EVIL")
        assert response.status_code == 400

    def test_rejects_unsupported_extension(self, client):
        response = self._upload(client, filename="payload.exe", genome="hg19")
        assert response.status_code == 400

    def test_rejects_filename_that_sanitizes_to_nothing(self, client):
        response = self._upload(client, filename="..", genome="hg19")
        assert response.status_code == 400

    def test_truncates_long_name(self, client, db):
        from mutagene.webapp.server import MAX_NAME_LENGTH

        response = self._upload(client, genome="hg19", name="A" * 5000)
        assert response.status_code == 200
        aid = response.get_json()["analysis_id"]
        assert len(db.get_analysis(aid)["name"]) == MAX_NAME_LENGTH

    def test_checkbox_on_enables_option(self, client, db):
        # Browsers submit "on" for a checked box with no value attribute.
        response = self._upload(client, genome="hg19", classify="on")
        aid = response.get_json()["analysis_id"]
        config = db.get_analysis(aid)["config"]
        assert config["classify"] is True

    def test_omitted_checkbox_disables_option(self, client, db):
        response = self._upload(client, genome="hg19")
        aid = response.get_json()["analysis_id"]
        config = db.get_analysis(aid)["config"]
        assert config["classify"] is False
        assert config["cluster"] is False


class TestRequestGuards:
    def test_rejects_foreign_host_header(self, client):
        response = client.get("/", headers={"Host": "evil.example.com"})
        assert response.status_code == 400

    def test_allows_loopback_host(self, client):
        assert client.get("/", headers={"Host": "127.0.0.1:5000"}).status_code == 200

    def test_blocks_cross_origin_state_change(self, client, db):
        aid = db.create_analysis("test.maf", "hg19")
        response = client.delete(
            f"/api/delete/{aid}", headers={"Origin": "http://evil.example.com"}
        )
        assert response.status_code == 403
        assert db.get_analysis(aid) is not None

    def test_allows_same_origin_state_change(self, client, db):
        aid = db.create_analysis("test.maf", "hg19")
        response = client.delete(f"/api/delete/{aid}", headers={"Origin": "http://localhost:5000"})
        assert response.status_code == 200

    def test_get_is_not_origin_checked(self, client):
        response = client.get("/", headers={"Origin": "http://evil.example.com"})
        assert response.status_code == 200


class TestDeleteSafety:
    def test_refuses_to_unlink_outside_managed_dirs(self, client, db, tmp_path):
        aid = db.create_analysis("test.maf", "hg19")
        outside = tmp_path / "outside" / "important.txt"
        outside.parent.mkdir(parents=True)
        outside.write_text("keep me")
        db.register_file(aid, "output", "important.txt", str(outside))

        assert client.delete(f"/api/delete/{aid}").status_code == 200
        assert outside.exists(), "file outside upload/results dirs must not be deleted"

    def test_deletes_managed_file(self, client, db, app):
        aid = db.create_analysis("test.maf", "hg19")
        managed = app.config["UPLOAD_FOLDER"] / "input.maf"
        managed.write_text("data")
        db.register_file(aid, "input_maf", "input.maf", str(managed))

        assert client.delete(f"/api/delete/{aid}").status_code == 200
        assert not managed.exists()

    def test_refuses_to_delete_running_analysis(self, client, db):
        aid = db.create_analysis("test.maf", "hg19")
        db.update_analysis_status(aid, "running")
        assert client.delete(f"/api/delete/{aid}").status_code == 409


class TestAnalyzeConcurrency:
    def test_second_start_is_rejected(self, client, db, app):
        from unittest.mock import MagicMock

        app.config["ANALYSIS_EXECUTOR"] = MagicMock()
        aid = db.create_analysis("test.maf", "hg19")

        assert client.post(f"/api/analyze/{aid}").status_code == 200
        assert client.post(f"/api/analyze/{aid}").status_code == 400
        # Only one worker may ever be submitted for a single claim.
        assert app.config["ANALYSIS_EXECUTOR"].submit.call_count == 1


class TestStartServerGuards:
    def test_debug_on_public_host_is_refused(self):
        from mutagene.webapp.server import start_server

        # The Werkzeug debugger is an RCE console; it must never be network-reachable.
        with pytest.raises(ValueError, match="debug mode"):
            start_server(host="0.0.0.0", port=5000, debug=True, open_browser=False)

    def test_debug_on_loopback_is_allowed_through_guard(self, monkeypatch):
        import mutagene.webapp.server as server_mod

        called = {}

        def fake_create_app(config=None):
            called["config"] = config
            raise RuntimeError("stop before binding")

        monkeypatch.setattr(server_mod, "create_app", fake_create_app)
        with pytest.raises(RuntimeError, match="stop before binding"):
            server_mod.start_server(host="127.0.0.1", port=8080, debug=True, open_browser=False)

        # Origins must track the actual port, not a hardcoded 5000.
        assert "http://127.0.0.1:8080" in called["config"]["ALLOWED_ORIGINS"]
