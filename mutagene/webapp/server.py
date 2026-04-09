"""Main Flask application for mutagene local web interface."""

import logging
import os
import sys
import threading
import webbrowser
from pathlib import Path

from flask import Flask, jsonify, render_template, request, send_file
from flask_socketio import SocketIO
from werkzeug.utils import secure_filename

from .database import DatabaseManager
from .genome_manager import GenomeManager

logger = logging.getLogger(__name__)


def create_app(config=None):
    """Create and configure the Flask application.

    Args:
        config: Optional configuration dictionary

    Returns:
        Flask application instance
    """
    # Create Flask app
    app = Flask(
        __name__,
        template_folder=str(Path(__file__).parent / "templates"),
        static_folder=str(Path(__file__).parent / "static"),
    )

    # Configuration
    import secrets

    app.config["SECRET_KEY"] = os.environ.get("SECRET_KEY") or secrets.token_hex(32)
    app.config["MAX_CONTENT_LENGTH"] = 500 * 1024 * 1024  # 500 MB max upload
    app.config["UPLOAD_FOLDER"] = Path.home() / ".mutagene" / "uploads"
    app.config["UPLOAD_FOLDER"].mkdir(parents=True, exist_ok=True)

    if config:
        app.config.update(config)

    # Initialize SocketIO for real-time updates
    socketio = SocketIO(
        app, cors_allowed_origins=["http://localhost:5000", "http://127.0.0.1:5000"]
    )

    # Initialize database
    db = DatabaseManager()

    # Initialize genome manager and check for required genomes
    genome_manager = GenomeManager()
    logger.info("Checking for required genome files...")

    # First check what's available
    available = genome_manager.get_available_genomes()
    if "hg19" in available and "hg38" in available:
        logger.info(f"All required genomes available: {', '.join(available)}")
    else:
        # Auto-download missing genomes
        logger.info("Downloading required genome files (this may take a few minutes)...")
        genome_status = genome_manager.check_and_download_required_genomes(
            required=["hg19", "hg38"], auto_download=True
        )

        if genome_status["downloaded"]:
            logger.info(f"Downloaded: {', '.join(genome_status['downloaded'])}")

        if genome_status["missing"]:
            logger.error(f"Failed to download genomes: {', '.join(genome_status['missing'])}")
            logger.error("Server may not function properly without genome files")
        else:
            logger.info("All required genomes are now available")

    # Store genome manager in app config for access in routes
    app.config["GENOME_MANAGER"] = genome_manager

    # Register routes
    register_routes(app, db, socketio)
    register_socketio_handlers(socketio, db)

    return app, socketio


def register_routes(app, db, socketio):
    """Register Flask routes."""

    @app.route("/")
    def index():
        """Home page with upload form."""
        recent_analyses = db.list_analyses(limit=10)
        return render_template("index.html", analyses=recent_analyses)

    @app.route("/history")
    def history():
        """Analysis history page."""
        analyses = db.list_analyses(limit=50)
        return render_template("history.html", analyses=analyses)

    @app.route("/api/upload", methods=["POST"])
    def upload_file():
        """Handle file upload."""
        if "file" not in request.files:
            return jsonify({"error": "No file provided"}), 400

        file = request.files["file"]
        if file.filename == "":
            return jsonify({"error": "Empty filename"}), 400

        # Get parameters
        name = request.form.get("name", file.filename)
        genome = request.form.get("genome", "hg19")
        signatures = request.form.get("signatures", "COSMICv3")

        # Save uploaded file with unique prefix to avoid collisions
        import uuid

        filename = f"{uuid.uuid4().hex[:8]}_{secure_filename(file.filename)}"
        upload_path = app.config["UPLOAD_FOLDER"] / filename
        file.save(str(upload_path))

        # Create analysis record
        analysis_id = db.create_analysis(
            name=name,
            genome=genome,
            signatures=signatures,
            config={
                "classify": request.form.get("classify", "true") == "true",
                "cluster": request.form.get("cluster", "true") == "true",
                "hotspots": request.form.get("hotspots", "true") == "true",
                "motifs": request.form.get("motifs", "true") == "true",
            },
        )

        # Register uploaded file
        db.register_file(analysis_id, "input_maf", filename, str(upload_path))

        return jsonify({"analysis_id": analysis_id, "message": "File uploaded successfully"})

    @app.route("/api/analyze/<int:analysis_id>", methods=["POST"])
    def start_analysis(analysis_id):
        """Start analysis in background thread."""
        analysis = db.get_analysis(analysis_id)
        if not analysis:
            return jsonify({"error": "Analysis not found"}), 404

        if analysis["status"] == "running":
            return jsonify({"error": "Analysis already running"}), 400

        # Update status
        db.update_analysis_status(analysis_id, "running")

        # Start analysis in background
        thread = threading.Thread(
            target=run_analysis, args=(analysis_id, db, socketio), daemon=True
        )
        thread.start()

        return jsonify({"status": "started", "analysis_id": analysis_id})

    @app.route("/api/analysis/<int:analysis_id>")
    def get_analysis(analysis_id):
        """Get analysis status and details."""
        analysis = db.get_analysis(analysis_id)
        if not analysis:
            return jsonify({"error": "Analysis not found"}), 404
        # Strip internal traceback from client response
        if analysis.get("error_message"):
            analysis["error_message"] = analysis["error_message"].split("\n")[0]
        return jsonify(analysis)

    @app.route("/results/<int:analysis_id>")
    def view_results(analysis_id):
        """View analysis results."""
        analysis = db.get_analysis(analysis_id)
        if not analysis:
            return "Analysis not found", 404

        if analysis["status"] != "complete":
            return render_template("analyzing.html", analysis=analysis)

        # Get all results
        results = db.get_all_results(analysis_id)
        files = db.get_all_files(analysis_id)

        return render_template("results.html", analysis=analysis, results=results, files=files)

    @app.route("/api/genomes")
    def get_genomes():
        """Get genome download status."""
        genome_manager = app.config["GENOME_MANAGER"]
        return jsonify(genome_manager.get_genome_info())

    @app.route("/api/genomes/<genome>/download", methods=["POST"])
    def download_genome(genome):
        """Download a genome reference file."""
        genome_manager = app.config["GENOME_MANAGER"]

        if genome not in genome_manager.SUPPORTED_GENOMES:
            return jsonify({"error": f"Unsupported genome: {genome}"}), 400

        if genome_manager.is_downloaded(genome):
            return jsonify({"message": f"{genome} already downloaded"})

        # Start download in background thread
        def download_task():
            try:
                socketio.emit("genome_download", {"genome": genome, "status": "downloading"})
                success = genome_manager.download_genome(genome)
                socketio.emit(
                    "genome_download",
                    {"genome": genome, "status": "complete" if success else "failed"},
                )
            except Exception as e:
                logger.error(f"Genome download failed: {e}")
                socketio.emit(
                    "genome_download", {"genome": genome, "status": "failed", "error": str(e)}
                )

        thread = threading.Thread(target=download_task, daemon=True)
        thread.start()

        return jsonify({"status": "started", "genome": genome})

    @app.route("/api/download/<int:file_id>")
    def download_file(file_id):
        """Download a file by ID."""
        file_record = db.get_file(file_id)
        if not file_record:
            return jsonify({"error": "File not found"}), 404

        file_path = Path(file_record["path"]).resolve()
        allowed_dirs = [
            Path(app.config["UPLOAD_FOLDER"]).resolve(),
            (Path.home() / ".mutagene" / "results").resolve(),
        ]
        if not any(str(file_path).startswith(str(d)) for d in allowed_dirs):
            return jsonify({"error": "Access denied"}), 403

        if not file_path.exists():
            return jsonify({"error": "File no longer exists on disk"}), 404

        return send_file(
            str(file_path),
            as_attachment=True,
            download_name=file_record["filename"],
        )

    @app.route("/api/delete/<int:analysis_id>", methods=["DELETE"])
    def delete_analysis(analysis_id):
        """Delete an analysis and associated files."""
        import shutil

        # Refuse to delete a running analysis
        analysis = db.get_analysis(analysis_id)
        if analysis and analysis.get("status") == "running":
            return jsonify({"error": "Cannot delete a running analysis"}), 409

        # Clean up files on disk
        files = db.get_all_files(analysis_id)
        for f in files:
            file_path = Path(f["path"])
            if file_path.exists():
                file_path.unlink()

        # Clean up result directory
        result_dir = Path.home() / ".mutagene" / "results" / str(analysis_id)
        if result_dir.exists():
            shutil.rmtree(result_dir)

        db.delete_analysis(analysis_id)
        return jsonify({"status": "deleted"})


def register_socketio_handlers(socketio, db):
    """Register WebSocket handlers."""

    @socketio.on("connect")
    def handle_connect():
        """Handle client connection."""
        logger.debug("Client connected")

    @socketio.on("disconnect")
    def handle_disconnect():
        """Handle client disconnection."""
        logger.debug("Client disconnected")


def run_analysis(analysis_id, db, socketio):
    """Run analysis in background thread.

    Args:
        analysis_id: Analysis ID
        db: Database manager
        socketio: SocketIO instance for progress updates
    """
    from pathlib import Path

    from .analysis import run_cohort_analysis

    try:
        # Get analysis details
        analysis = db.get_analysis(analysis_id)
        if not analysis:
            raise ValueError(f"Analysis {analysis_id} not found")

        # Get input file
        files = db.get_all_files(analysis_id)
        input_files = [f for f in files if f["file_type"] == "input_maf"]
        if not input_files:
            raise ValueError("No input file found")

        input_path = Path(input_files[0]["path"])
        output_dir = Path.home() / ".mutagene" / "results" / str(analysis_id)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Emit progress updates
        socketio.emit(
            "progress", {"analysis_id": analysis_id, "step": "Reading input file", "progress": 10}
        )

        # Run analysis
        results = run_cohort_analysis(
            input_file=input_path,
            output_dir=output_dir,
            genome=analysis["genome"],
            signatures_set=analysis.get("signatures", "COSMICv3"),
            config=analysis.get("config", {}),
        )

        # Update analysis with results
        with db.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                "UPDATE analyses SET samples = ?, mutations = ? WHERE id = ?",
                (results["samples"], results["mutations"], analysis_id),
            )
            conn.commit()

        # Register output files
        for file_info in results["files"]:
            db.register_file(
                analysis_id, file_info["type"], Path(file_info["path"]).name, file_info["path"]
            )

        # Store analysis results in database
        if results.get("profiles"):
            db.store_result(analysis_id, "profile", results["profiles"])

        if results.get("signatures"):
            db.store_result(analysis_id, "signatures", results["signatures"])

        if results.get("classification"):
            db.store_result(analysis_id, "classification", results["classification"])

        if results.get("drivers"):
            db.store_result(analysis_id, "drivers", results["drivers"])

        if results.get("clustering"):
            db.store_result(analysis_id, "clustering", results["clustering"])

        if results.get("genome_warning"):
            db.store_result(analysis_id, "genome_warning", results["genome_warning"])

        # Mark as complete
        db.update_analysis_status(analysis_id, "complete")
        socketio.emit("complete", {"analysis_id": analysis_id})

    except Exception as e:
        import traceback

        error_msg = f"{str(e)}\n{traceback.format_exc()}"
        db.update_analysis_status(analysis_id, "error", error_msg)
        socketio.emit("error", {"analysis_id": analysis_id, "error": str(e)})


def start_server(host="127.0.0.1", port=5000, debug=False, open_browser=True):
    """Start the Flask development server.

    Args:
        host: Host to bind to
        port: Port to bind to
        debug: Enable debug mode
        open_browser: Automatically open browser
    """
    app, socketio = create_app()

    print(
        f"""
╔══════════════════════════════════════════════════════════════╗
║                    MutaGene Local Server                     ║
╠══════════════════════════════════════════════════════════════╣
║ Analysis engine initialized                                  ║
║ Database ready: ~/.mutagene/results.db                       ║
║ Server running at http://{host}:{port:<5}                       ║
║                                                              ║
║ Open your browser to get started!                           ║
║ Press Ctrl+C to stop                                         ║
╚══════════════════════════════════════════════════════════════╝
"""
    )

    if open_browser:
        # Open browser after short delay
        def open_browser_delayed():
            import time

            time.sleep(1.5)
            webbrowser.open(f"http://{host}:{port}")

        threading.Thread(target=open_browser_delayed, daemon=True).start()

    # Run with use_reloader=False and let Werkzeug handle signals naturally
    try:
        socketio.run(
            app,
            host=host,
            port=port,
            debug=debug,
            use_reloader=False,
            log_output=not debug,  # Reduce log spam in production
        )
    except KeyboardInterrupt:
        print("\n\nShutting down server...")
        sys.exit(0)
