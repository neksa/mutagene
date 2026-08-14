"""Re-running an analysis must not damage the previous run's output.

A re-run writes into a staging directory and swaps it into place only after the
analysis succeeds and its results commit, so a failure leaves both the files on
disk and the rows in the database as they were.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from mutagene.webapp.database.manager import DatabaseManager
from mutagene.webapp.server import _recover_interrupted_runs, run_analysis


@pytest.fixture
def db(tmp_path):
    return DatabaseManager(db_path=tmp_path / "test.db")


def _seed_analysis(db, tmp_path):
    """An analysis with an uploaded input and one completed run already stored."""
    analysis_id = db.create_analysis("s.maf", "hg38")
    maf = tmp_path / "s.maf"
    maf.write_text("input")
    db.register_file(analysis_id, "input_maf", "s.maf", str(maf))

    results_folder = tmp_path / "results"
    previous = results_folder / str(analysis_id)
    previous.mkdir(parents=True)
    (previous / "profile.tsv").write_text("PREVIOUS")

    db.register_file(analysis_id, "profile", "profile.tsv", str(previous / "profile.tsv"))
    db.store_result(analysis_id, "profile", {"run": "previous"})
    return analysis_id, results_folder


def _fake_run(payload):
    """Stand in for run_cohort_analysis, writing payload into the given output_dir."""

    def _run(input_file, output_dir, genome, signatures_set, config):
        out = output_dir / "profile.tsv"
        out.write_text(payload)
        return {
            "samples": 1,
            "mutations": 5,
            "files": [{"type": "profile", "path": str(out), "size": out.stat().st_size}],
            "profiles": {"run": payload},
        }

    return _run


def test_successful_rerun_replaces_output(db, tmp_path):
    analysis_id, results_folder = _seed_analysis(db, tmp_path)

    with patch("mutagene.webapp.analysis.run_cohort_analysis", _fake_run("NEW")):
        run_analysis(analysis_id, db, MagicMock(), results_folder)

    final_dir = results_folder / str(analysis_id)
    assert (final_dir / "profile.tsv").read_text() == "NEW"

    rows = db.get_results_by_type(analysis_id, "profile")
    assert len(rows) == 1, "the stale row must not survive alongside the new one"
    assert rows[0]["data"] == {"run": "NEW"}
    assert db.get_analysis(analysis_id)["status"] == "complete"

    # No staging or backup directories left behind.
    leftovers = [p.name for p in results_folder.iterdir() if p.name.startswith(".")]
    assert leftovers == []


def test_failed_rerun_leaves_previous_output_intact(db, tmp_path):
    analysis_id, results_folder = _seed_analysis(db, tmp_path)

    def _boom(input_file, output_dir, genome, signatures_set, config):
        # Write output first, then fail: the half-written run must not escape staging.
        (output_dir / "profile.tsv").write_text("HALF-WRITTEN")
        raise RuntimeError("Variant allele is not defined in MAF file")

    with patch("mutagene.webapp.analysis.run_cohort_analysis", _boom):
        run_analysis(analysis_id, db, MagicMock(), results_folder)

    final_dir = results_folder / str(analysis_id)
    assert (
        final_dir / "profile.tsv"
    ).read_text() == "PREVIOUS", "a failed re-run overwrote the previous run's file on disk"

    rows = db.get_results_by_type(analysis_id, "profile")
    assert len(rows) == 1
    assert rows[0]["data"] == {"run": "previous"}
    assert db.get_analysis(analysis_id)["status"] == "error"

    leftovers = [p.name for p in results_folder.iterdir() if p.name.startswith(".")]
    assert leftovers == [], "staging directory was not cleaned up"


def test_first_run_needs_no_previous_directory(db, tmp_path):
    analysis_id = db.create_analysis("s.maf", "hg38")
    maf = tmp_path / "s.maf"
    maf.write_text("input")
    db.register_file(analysis_id, "input_maf", "s.maf", str(maf))
    results_folder = tmp_path / "results"

    with patch("mutagene.webapp.analysis.run_cohort_analysis", _fake_run("FIRST")):
        run_analysis(analysis_id, db, MagicMock(), results_folder)

    assert (results_folder / str(analysis_id) / "profile.tsv").read_text() == "FIRST"
    assert db.get_analysis(analysis_id)["status"] == "complete"


def test_malformed_results_leave_previous_output_untouched(db, tmp_path):
    """A bad results payload must fail before any file is moved."""
    analysis_id, results_folder = _seed_analysis(db, tmp_path)

    def _malformed(input_file, output_dir, genome, signatures_set, config):
        (output_dir / "profile.tsv").write_text("NEW")
        return {"samples": 1, "mutations": 5}  # missing "files"

    with patch("mutagene.webapp.analysis.run_cohort_analysis", _malformed):
        run_analysis(analysis_id, db, MagicMock(), results_folder)

    final_dir = results_folder / str(analysis_id)
    assert (final_dir / "profile.tsv").read_text() == "PREVIOUS"
    assert db.get_results_by_type(analysis_id, "profile")[0]["data"] == {"run": "previous"}
    assert [p.name for p in results_folder.iterdir() if p.name.startswith(".")] == []


def test_failed_publish_restores_previous_files_and_rows(db, tmp_path):
    """If the database write fails, the swapped-in files are rolled back too."""
    analysis_id, results_folder = _seed_analysis(db, tmp_path)

    with (
        patch("mutagene.webapp.analysis.run_cohort_analysis", _fake_run("NEW")),
        patch.object(
            DatabaseManager, "publish_run_results", side_effect=RuntimeError("db exploded")
        ),
    ):
        run_analysis(analysis_id, db, MagicMock(), results_folder)

    final_dir = results_folder / str(analysis_id)
    assert (
        final_dir / "profile.tsv"
    ).read_text() == "PREVIOUS", "the previous run's files must come back when the publish fails"
    assert db.get_results_by_type(analysis_id, "profile")[0]["data"] == {"run": "previous"}
    assert db.get_analysis(analysis_id)["status"] == "error"
    assert [p.name for p in results_folder.iterdir() if p.name.startswith(".")] == []


def test_nested_output_paths_keep_their_subdirectory(db, tmp_path):
    analysis_id, results_folder = _seed_analysis(db, tmp_path)

    def _nested(input_file, output_dir, genome, signatures_set, config):
        sub = output_dir / "plots"
        sub.mkdir()
        out = sub / "sig.png"
        out.write_text("img")
        return {
            "samples": 1,
            "mutations": 5,
            "files": [{"type": "plot", "path": str(out), "size": 3}],
            "profiles": {"run": "nested"},
        }

    with patch("mutagene.webapp.analysis.run_cohort_analysis", _nested):
        run_analysis(analysis_id, db, MagicMock(), results_folder)

    registered = db.get_files_by_type(analysis_id, "plot")
    assert len(registered) == 1
    registered_path = Path(registered[0]["path"])
    assert registered_path.exists(), "registered path must point at the real file"
    assert registered_path.parent.name == "plots"


class TestCrashRecovery:
    """A crash between the directory swap and the commit must be reconciled."""

    def test_discards_staging_left_by_a_crash(self, db, tmp_path):
        results_folder = tmp_path / "results"
        (results_folder / ".7.incoming").mkdir(parents=True)
        (results_folder / ".7.incoming" / "half.tsv").write_text("half")

        restored, discarded = _recover_interrupted_runs(results_folder, db)

        assert (restored, discarded) == (0, 1)
        assert not (results_folder / ".7.incoming").exists()

    def test_restores_backup_when_the_commit_did_not_land(self, db, tmp_path):
        analysis_id = db.create_analysis("s.maf", "hg38")
        db.update_analysis_status(analysis_id, "error", "interrupted")

        results_folder = tmp_path / "results"
        final_dir = results_folder / str(analysis_id)
        final_dir.mkdir(parents=True)
        (final_dir / "profile.tsv").write_text("UNCOMMITTED NEW")
        backup = results_folder / f".{analysis_id}.previous"
        backup.mkdir()
        (backup / "profile.tsv").write_text("PREVIOUS")

        restored, discarded = _recover_interrupted_runs(results_folder, db)

        assert (restored, discarded) == (1, 0)
        assert (final_dir / "profile.tsv").read_text() == "PREVIOUS"
        assert not backup.exists()

    def test_keeps_new_output_when_the_commit_landed(self, db, tmp_path):
        # Crash after the commit but before backup cleanup: restoring here would
        # silently revert a run that actually succeeded.
        analysis_id = db.create_analysis("s.maf", "hg38")
        db.update_analysis_status(analysis_id, "complete")

        results_folder = tmp_path / "results"
        final_dir = results_folder / str(analysis_id)
        final_dir.mkdir(parents=True)
        (final_dir / "profile.tsv").write_text("COMMITTED NEW")
        backup = results_folder / f".{analysis_id}.previous"
        backup.mkdir()
        (backup / "profile.tsv").write_text("PREVIOUS")

        restored, discarded = _recover_interrupted_runs(results_folder, db)

        assert (restored, discarded) == (0, 1)
        assert (final_dir / "profile.tsv").read_text() == "COMMITTED NEW"
        assert not backup.exists()

    def test_ignores_unrelated_hidden_entries(self, db, tmp_path):
        results_folder = tmp_path / "results"
        results_folder.mkdir(parents=True)
        (results_folder / ".DS_Store").write_text("junk")
        (results_folder / ".notanid.previous").mkdir()

        restored, discarded = _recover_interrupted_runs(results_folder, db)

        assert (results_folder / ".DS_Store").exists()
        # A non-numeric name has no analysis to consult; it is restored, not lost.
        assert (restored, discarded) == (1, 0)


def test_reporting_failure_does_not_undo_a_committed_run(db, tmp_path):
    """A crash after the commit must not relabel a successful run as failed."""
    analysis_id, results_folder = _seed_analysis(db, tmp_path)

    def _emit(event, *args, **kwargs):
        # Only the post-commit notification fails; progress updates succeed.
        if event == "complete":
            raise RuntimeError("socket gone")

    socketio = MagicMock()
    socketio.emit.side_effect = _emit

    with patch("mutagene.webapp.analysis.run_cohort_analysis", _fake_run("NEW")):
        run_analysis(analysis_id, db, socketio, results_folder)

    assert (
        db.get_analysis(analysis_id)["status"] == "complete"
    ), "the committed run must stay complete when only reporting failed"
    assert db.get_results_by_type(analysis_id, "profile")[0]["data"] == {"run": "NEW"}
    assert (results_folder / str(analysis_id) / "profile.tsv").read_text() == "NEW"

    # Recovery must agree with the committed rows rather than roll them back.
    _recover_interrupted_runs(results_folder, db)
    assert (results_folder / str(analysis_id) / "profile.tsv").read_text() == "NEW"
