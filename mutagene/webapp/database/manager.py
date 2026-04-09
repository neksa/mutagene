"""Database manager for mutagene webapp."""

import sqlite3
from contextlib import contextmanager
from pathlib import Path

from .models import Analysis, File, Result, init_db


class DatabaseManager:
    """Manages database connections and operations."""

    def __init__(self, db_path=None):
        """Initialize database manager.

        Args:
            db_path: Path to SQLite database. Defaults to ~/.mutagene/results.db
        """
        if db_path is None:
            db_path = Path.home() / ".mutagene" / "results.db"
        else:
            db_path = Path(db_path)

        # Create directory if it doesn't exist
        db_path.parent.mkdir(parents=True, exist_ok=True)

        self.db_path = str(db_path)

        # Initialize schema if database doesn't exist
        if not db_path.exists():
            init_db(self.db_path)

    @contextmanager
    def get_connection(self):
        """Context manager for database connections."""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        try:
            yield conn
        finally:
            conn.close()

    # Analysis methods
    def create_analysis(self, name, genome, signatures="COSMICv3", config=None):
        """Create new analysis."""
        with self.get_connection() as conn:
            return Analysis.create(conn, name, genome, signatures, config)

    def get_analysis(self, analysis_id):
        """Get analysis by ID."""
        with self.get_connection() as conn:
            return Analysis.get(conn, analysis_id)

    def list_analyses(self, limit=50):
        """List all analyses."""
        with self.get_connection() as conn:
            return Analysis.list_all(conn, limit)

    def update_analysis_status(self, analysis_id, status, error_message=None):
        """Update analysis status."""
        with self.get_connection() as conn:
            Analysis.update_status(conn, analysis_id, status, error_message)

    def update_analysis_counts(self, analysis_id, samples, mutations):
        """Update sample and mutation counts."""
        with self.get_connection() as conn:
            Analysis.update_counts(conn, analysis_id, samples, mutations)

    def delete_analysis(self, analysis_id):
        """Delete analysis."""
        with self.get_connection() as conn:
            Analysis.delete(conn, analysis_id)

    # Result methods
    def store_result(self, analysis_id, result_type, data, sample_id=None):
        """Store analysis results."""
        with self.get_connection() as conn:
            return Result.create(conn, analysis_id, result_type, data, sample_id)

    def get_results_by_type(self, analysis_id, result_type):
        """Get results of specific type."""
        with self.get_connection() as conn:
            return Result.get_by_type(conn, analysis_id, result_type)

    def get_all_results(self, analysis_id):
        """Get all results for analysis."""
        with self.get_connection() as conn:
            return Result.get_all(conn, analysis_id)

    # File methods
    def register_file(self, analysis_id, file_type, filename, path):
        """Register a file."""
        with self.get_connection() as conn:
            return File.create(conn, analysis_id, file_type, filename, path)

    def get_file(self, file_id):
        """Get a file by ID."""
        with self.get_connection() as conn:
            return File.get_by_id(conn, file_id)

    def get_files_by_type(self, analysis_id, file_type):
        """Get files of specific type."""
        with self.get_connection() as conn:
            return File.get_by_type(conn, analysis_id, file_type)

    def get_all_files(self, analysis_id):
        """Get all files for analysis."""
        with self.get_connection() as conn:
            return File.get_all(conn, analysis_id)
