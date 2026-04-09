"""SQLite database models for storing analysis results."""

import json
import sqlite3
from pathlib import Path


class _SafeEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy types gracefully."""

    def default(self, obj):
        try:
            import numpy as np

            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
        except ImportError:
            pass
        return super().default(obj)


def init_db(db_path):
    """Initialize the database schema.

    Args:
        db_path: Path to SQLite database file
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Analyses table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS analyses (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL,
            uploaded_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            genome TEXT NOT NULL,
            signatures TEXT,
            status TEXT DEFAULT 'pending',
            samples INTEGER DEFAULT 0,
            mutations INTEGER DEFAULT 0,
            error_message TEXT,
            config JSON
        )
    """)

    # Results table - stores analysis results as JSON
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            analysis_id INTEGER NOT NULL,
            result_type TEXT NOT NULL,
            sample_id TEXT,
            data JSON NOT NULL,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (analysis_id) REFERENCES analyses(id) ON DELETE CASCADE
        )
    """)

    # Files table - tracks input/output files
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS files (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            analysis_id INTEGER NOT NULL,
            file_type TEXT NOT NULL,
            filename TEXT NOT NULL,
            path TEXT NOT NULL,
            size_bytes INTEGER,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (analysis_id) REFERENCES analyses(id) ON DELETE CASCADE
        )
    """)

    # Create indices
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_analyses_status ON analyses(status)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_results_analysis ON results(analysis_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_results_type ON results(result_type)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_files_analysis ON files(analysis_id)")

    conn.commit()
    conn.close()


class Analysis:
    """Represents an analysis run."""

    @staticmethod
    def create(conn, name, genome, signatures="COSMICv3", config=None):
        """Create a new analysis record.

        Args:
            conn: SQLite connection
            name: Analysis name
            genome: Genome build (hg19, hg38, etc.)
            signatures: Signature set to use
            config: Configuration dictionary

        Returns:
            analysis_id: ID of created analysis
        """
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO analyses (name, genome, signatures, config)
            VALUES (?, ?, ?, ?)
        """,
            (name, genome, signatures, json.dumps(config) if config else None),
        )
        conn.commit()
        return cursor.lastrowid

    @staticmethod
    def get(conn, analysis_id):
        """Get analysis by ID."""
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM analyses WHERE id = ?", (analysis_id,))
        row = cursor.fetchone()
        if row:
            result = dict(zip([col[0] for col in cursor.description], row))
            # Deserialize JSON config field if present
            if result.get("config"):
                result["config"] = json.loads(result["config"])
            return result
        return None

    @staticmethod
    def update_status(conn, analysis_id, status, error_message=None):
        """Update analysis status."""
        cursor = conn.cursor()
        cursor.execute(
            """
            UPDATE analyses
            SET status = ?, error_message = ?
            WHERE id = ?
        """,
            (status, error_message, analysis_id),
        )
        conn.commit()

    @staticmethod
    def update_counts(conn, analysis_id, samples, mutations):
        """Update sample and mutation counts."""
        cursor = conn.cursor()
        cursor.execute(
            """
            UPDATE analyses
            SET samples = ?, mutations = ?
            WHERE id = ?
        """,
            (samples, mutations, analysis_id),
        )
        conn.commit()

    @staticmethod
    def list_all(conn, limit=50):
        """List all analyses, most recent first."""
        cursor = conn.cursor()
        cursor.execute(
            """
            SELECT * FROM analyses
            ORDER BY uploaded_at DESC
            LIMIT ?
        """,
            (limit,),
        )
        columns = [col[0] for col in cursor.description]
        results = []
        for row in cursor.fetchall():
            result = dict(zip(columns, row))
            # Deserialize JSON config field if present
            if result.get("config"):
                result["config"] = json.loads(result["config"])
            results.append(result)
        return results

    @staticmethod
    def delete(conn, analysis_id):
        """Delete an analysis and all related data."""
        cursor = conn.cursor()
        cursor.execute("DELETE FROM analyses WHERE id = ?", (analysis_id,))
        conn.commit()


class Result:
    """Represents analysis results."""

    @staticmethod
    def create(conn, analysis_id, result_type, data, sample_id=None):
        """Store analysis results.

        Args:
            conn: SQLite connection
            analysis_id: Analysis ID
            result_type: Type of result (profile, signature, classification, etc.)
            data: Result data (will be JSON serialized)
            sample_id: Optional sample identifier
        """
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO results (analysis_id, result_type, sample_id, data)
            VALUES (?, ?, ?, ?)
        """,
            (analysis_id, result_type, sample_id, json.dumps(data, cls=_SafeEncoder)),
        )
        conn.commit()
        return cursor.lastrowid

    @staticmethod
    def get_by_type(conn, analysis_id, result_type):
        """Get all results of a specific type for an analysis."""
        cursor = conn.cursor()
        cursor.execute(
            """
            SELECT * FROM results
            WHERE analysis_id = ? AND result_type = ?
            ORDER BY created_at
        """,
            (analysis_id, result_type),
        )
        columns = [col[0] for col in cursor.description]
        results = []
        for row in cursor.fetchall():
            result = dict(zip(columns, row))
            result["data"] = json.loads(result["data"])
            results.append(result)
        return results

    @staticmethod
    def get_all(conn, analysis_id):
        """Get all results for an analysis."""
        cursor = conn.cursor()
        cursor.execute(
            """
            SELECT * FROM results
            WHERE analysis_id = ?
            ORDER BY result_type, created_at
        """,
            (analysis_id,),
        )
        columns = [col[0] for col in cursor.description]
        results = []
        for row in cursor.fetchall():
            result = dict(zip(columns, row))
            result["data"] = json.loads(result["data"])
            results.append(result)
        return results


class File:
    """Represents files associated with analyses."""

    @staticmethod
    def create(conn, analysis_id, file_type, filename, path):
        """Register a file.

        Args:
            conn: SQLite connection
            analysis_id: Analysis ID
            file_type: Type of file (input_maf, profile_tsv, plot_png, etc.)
            filename: Original filename
            path: Full path to file
        """
        size_bytes = Path(path).stat().st_size if Path(path).exists() else 0
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO files (analysis_id, file_type, filename, path, size_bytes)
            VALUES (?, ?, ?, ?, ?)
        """,
            (analysis_id, file_type, filename, path, size_bytes),
        )
        conn.commit()
        return cursor.lastrowid

    @staticmethod
    def get_by_id(conn, file_id):
        """Get a file by ID."""
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM files WHERE id = ?", (file_id,))
        row = cursor.fetchone()
        if row:
            return dict(zip([col[0] for col in cursor.description], row))
        return None

    @staticmethod
    def get_by_type(conn, analysis_id, file_type):
        """Get files of a specific type."""
        cursor = conn.cursor()
        cursor.execute(
            """
            SELECT * FROM files
            WHERE analysis_id = ? AND file_type = ?
            ORDER BY created_at
        """,
            (analysis_id, file_type),
        )
        columns = [col[0] for col in cursor.description]
        return [dict(zip(columns, row)) for row in cursor.fetchall()]

    @staticmethod
    def get_all(conn, analysis_id):
        """Get all files for an analysis."""
        cursor = conn.cursor()
        cursor.execute(
            """
            SELECT * FROM files
            WHERE analysis_id = ?
            ORDER BY file_type, created_at
        """,
            (analysis_id,),
        )
        columns = [col[0] for col in cursor.description]
        return [dict(zip(columns, row)) for row in cursor.fetchall()]
