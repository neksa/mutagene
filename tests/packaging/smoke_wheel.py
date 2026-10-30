"""Run this script outside the checkout after installing the wheel."""

import subprocess
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

from mutagene.webapp.database import DatabaseManager
from mutagene.webapp.genome_manager import GenomeManager
from mutagene.webapp.server import create_app

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    subprocess.run(
        [sys.executable, "-m", "mutagene", "demo", "--output-dir", str(root / "demo")],
        check=True,
        cwd=root,
    )
    with (
        patch(
            "mutagene.webapp.server.DatabaseManager", return_value=DatabaseManager(root / "app.db")
        ),
        patch("mutagene.webapp.server.GenomeManager", return_value=GenomeManager(root / "genomes")),
    ):
        app, _ = create_app(
            {
                "TESTING": True,
                "AUTO_DOWNLOAD_GENOMES": False,
                "UPLOAD_FOLDER": root / "uploads",
                "RESULTS_FOLDER": root / "results",
            }
        )
        try:
            client = app.test_client()
            assert client.get("/").status_code == 200
            assert client.get("/static/img/mutagene_logo.png").status_code == 200
        finally:
            app.config["ANALYSIS_EXECUTOR"].shutdown()
print("Installed wheel: offline demo, homepage, and static asset pass")
