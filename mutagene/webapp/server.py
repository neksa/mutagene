"""Main Flask application for mutagene local web interface."""

import logging
import os
import shutil
import sys
import threading
import webbrowser
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from flask import Flask, jsonify, render_template, request, send_file
from flask_socketio import SocketIO
from werkzeug.utils import secure_filename

from .database import DatabaseManager
from .genome_manager import GenomeManager

logger = logging.getLogger(__name__)

# Snapshot the genome list at import so request validation is independent of the
# GenomeManager instance held in app config (which tests replace with a mock).
SUPPORTED_GENOMES = frozenset(GenomeManager.SUPPORTED_GENOMES)

# Signature sets accepted by mutagene.io.profile.read_signatures
SUPPORTED_SIGNATURE_SETS = frozenset({"MGA", "MGB", "COSMICv2", "COSMICv3", "KUCAB"})

# Mutation file formats the analysis pipeline can read, plus the archive/compression
# wrappers extract_input_file() understands.
ALLOWED_EXTENSIONS = frozenset({".maf", ".vcf", ".txt", ".tsv", ".gz", ".tgz"})

MAX_NAME_LENGTH = 200


def _has_allowed_extension(filename: str) -> bool:
    """True if *filename* ends in an accepted mutation-file extension."""
    return Path(filename).suffix.lower() in ALLOWED_EXTENSIONS


def _form_flag(field: str) -> bool:
    """Read a checkbox-style form field.

    An unchecked HTML checkbox is omitted from the submission entirely, and a
    checked one without an explicit value attribute submits the string "on" —
    so comparing against "true" with a "true" default inverts every option.
    """
    value = request.form.get(field)
    if value is None:
        return False
    return value.strip().lower() in ("true", "on", "1", "yes")


def _managed_dirs(app) -> list[Path]:
    """Directories this app is allowed to read from and delete within."""
    return [Path(app.config["UPLOAD_FOLDER"]), Path(app.config["RESULTS_FOLDER"])]


def _is_within(path: Path, directory: Path) -> bool:
    """True if *path* resolves to a location inside *directory*.

    Uses path-component comparison rather than string prefixing, which would
    treat "/data/uploads-evil" as living under "/data/uploads".
    """
    try:
        return path.resolve().is_relative_to(directory.resolve())
    except (OSError, ValueError):
        return False


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

    secret = os.environ.get("SECRET_KEY")
    if not secret:
        secret = secrets.token_hex(32)
        logger.info(
            "No SECRET_KEY set; using auto-generated key (sessions won't persist across restarts)"
        )
    app.config["SECRET_KEY"] = secret
    app.config["MAX_CONTENT_LENGTH"] = 500 * 1024 * 1024  # 500 MB max upload
    app.config["UPLOAD_FOLDER"] = Path.home() / ".mutagene" / "uploads"
    app.config["RESULTS_FOLDER"] = Path.home() / ".mutagene" / "results"
    # Origins allowed to make state-changing requests. start_server() extends this
    # with the actual bind host/port; the defaults cover the common case.
    app.config["ALLOWED_ORIGINS"] = [
        "http://localhost:5000",
        "http://127.0.0.1:5000",
    ]
    # Hostnames this server will answer to. None disables the check, which is
    # what start_server() does for a deliberate non-loopback bind where the
    # hostname clients use cannot be known ahead of time.
    app.config["ALLOWED_HOSTS"] = ["localhost", "127.0.0.1", "::1"]
    app.config["AUTO_DOWNLOAD_GENOMES"] = True
    app.config["MAX_CONCURRENT_ANALYSES"] = 2

    # Apply caller overrides before creating directories, so a test or embedder
    # pointing UPLOAD_FOLDER elsewhere does not also create the default location.
    if config:
        app.config.update(config)

    app.config["UPLOAD_FOLDER"] = Path(app.config["UPLOAD_FOLDER"])
    app.config["RESULTS_FOLDER"] = Path(app.config["RESULTS_FOLDER"])
    app.config["UPLOAD_FOLDER"].mkdir(parents=True, exist_ok=True)
    app.config["RESULTS_FOLDER"].mkdir(parents=True, exist_ok=True)

    # Initialize SocketIO for real-time updates
    socketio = SocketIO(app, cors_allowed_origins=list(app.config["ALLOWED_ORIGINS"]))

    # Initialize database
    db = DatabaseManager()

    # Analyses left mid-flight by a previous process can never finish; without
    # this they stay 'running' forever and can never be restarted or deleted.
    reclaimed = db.reset_stale_running()
    if reclaimed:
        logger.warning(f"Marked {reclaimed} interrupted analysis(es) as failed")

    # A crash between swapping a run's output into place and committing its rows
    # leaves the directory and the database disagreeing; put them back in step.
    restored, discarded = _recover_interrupted_runs(app.config["RESULTS_FOLDER"], db)
    if restored or discarded:
        logger.warning(
            f"Recovered interrupted analysis output: restored {restored} previous "
            f"result directory(ies), discarded {discarded} incomplete one(s)"
        )

    app.config["ANALYSIS_EXECUTOR"] = ThreadPoolExecutor(
        max_workers=app.config["MAX_CONCURRENT_ANALYSES"],
        thread_name_prefix="mutagene-analysis",
    )

    # Initialize genome manager and check for required genomes
    genome_manager = GenomeManager()
    app.config["GENOME_MANAGER"] = genome_manager

    _ensure_genomes(genome_manager, auto_download=app.config["AUTO_DOWNLOAD_GENOMES"])

    # Register routes
    register_request_guards(app)
    register_routes(app, db, socketio)
    register_socketio_handlers(socketio, db)

    return app, socketio


def _ensure_genomes(genome_manager, auto_download=True, required=("hg19", "hg38")):
    """Check for required genomes, downloading missing ones in the background.

    The reference genomes are several GB each. Downloading them inline would keep
    the server from accepting connections for many minutes with no UI feedback,
    so missing genomes are fetched on a daemon thread while the app starts.
    """
    logger.info("Checking for required genome files...")
    available = genome_manager.get_available_genomes()
    missing = [g for g in required if g not in available]

    if not missing:
        logger.info(f"All required genomes available: {', '.join(available)}")
        return

    if not auto_download:
        logger.warning(f"Missing genome files: {', '.join(missing)}")
        logger.warning("Download them with: mutagene fetch genome <name>")
        return

    def download_missing():
        logger.info(f"Downloading genome files in background: {', '.join(missing)}")
        status = genome_manager.check_and_download_required_genomes(
            required=list(missing), auto_download=True
        )
        if status["downloaded"]:
            logger.info(f"Downloaded: {', '.join(status['downloaded'])}")
        still_missing = [g for g in status["missing"] if g not in status["downloaded"]]
        if still_missing:
            logger.error(f"Failed to download genomes: {', '.join(still_missing)}")
            logger.error("Analyses requiring these assemblies will fail until they are present")

    threading.Thread(target=download_missing, daemon=True).start()


def register_request_guards(app):
    """Reject cross-origin and rebound-DNS requests to this local-only server.

    The app has no authentication and binds to localhost, so any page in the
    user's browser could otherwise drive it: a plain HTML form can POST
    multipart data to /api/upload cross-site, and a hostname resolving to
    127.0.0.1 could reach every endpoint. Checking Origin and Host closes both
    without adding a token-based CSRF dependency.
    """

    @app.before_request
    def guard_request():
        allowed = set(app.config.get("ALLOWED_ORIGINS", []))
        allowed_hosts = app.config.get("ALLOWED_HOSTS")

        # DNS-rebinding guard: only answer to hostnames we expect. Disabled
        # (ALLOWED_HOSTS = None) for a deliberate non-loopback bind.
        if allowed_hosts is not None:
            host = (request.host or "").rsplit(":", 1)[0].strip("[]")
            if host and host not in allowed_hosts:
                return jsonify({"error": "Invalid Host header"}), 400

        # CSRF guard: browsers always send Origin on cross-site state changes.
        if request.method not in ("GET", "HEAD", "OPTIONS"):
            origin = request.headers.get("Origin")
            if origin is not None and origin not in allowed:
                return jsonify({"error": "Cross-origin request blocked"}), 403
        return None

    @app.errorhandler(413)
    def handle_too_large(_error):
        limit_mb = app.config["MAX_CONTENT_LENGTH"] // (1024 * 1024)
        return jsonify({"error": f"File too large (limit {limit_mb} MB)"}), 413


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
        if not file.filename:
            return jsonify({"error": "Empty filename"}), 400

        # Validate parameters against known-good values before they reach the
        # filesystem (genome becomes part of a path) or the analysis pipeline.
        genome = request.form.get("genome", "hg19")
        if genome not in SUPPORTED_GENOMES:
            return jsonify({"error": f"Unsupported genome: {genome}"}), 400

        signatures = request.form.get("signatures", "COSMICv3")
        if signatures not in SUPPORTED_SIGNATURE_SETS:
            return jsonify({"error": f"Unsupported signature set: {signatures}"}), 400

        name = (request.form.get("name") or file.filename).strip()[:MAX_NAME_LENGTH]
        if not name:
            return jsonify({"error": "Empty analysis name"}), 400

        safe_name = secure_filename(file.filename)
        if not safe_name:
            return jsonify({"error": "Invalid filename"}), 400
        if not _has_allowed_extension(safe_name):
            return (
                jsonify(
                    {
                        "error": f"Unsupported file type. Allowed: {', '.join(sorted(ALLOWED_EXTENSIONS))}"
                    }
                ),
                400,
            )

        # Save uploaded file with unique prefix to avoid collisions
        import uuid

        filename = f"{uuid.uuid4().hex[:8]}_{safe_name}"
        upload_path = app.config["UPLOAD_FOLDER"] / filename
        file.save(str(upload_path))

        # Create analysis record
        analysis_id = db.create_analysis(
            name=name,
            genome=genome,
            signatures=signatures,
            config={
                "classify": _form_flag("classify"),
                "cluster": _form_flag("cluster"),
                "hotspots": _form_flag("hotspots"),
                "motifs": _form_flag("motifs"),
                "keep_filtered": _form_flag("keep_filtered"),
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

        # Claim the analysis atomically. Checking the status and then updating it
        # in two statements lets two concurrent POSTs both pass the check and
        # start duplicate threads writing the same output directory.
        if not db.try_claim_analysis(analysis_id):
            return jsonify({"error": "Analysis already running"}), 400

        # Run on a bounded pool: analyses are CPU and memory heavy, and an
        # unbounded thread per request lets a burst of uploads exhaust the box.
        app.config["ANALYSIS_EXECUTOR"].submit(
            run_analysis, analysis_id, db, socketio, app.config["RESULTS_FOLDER"]
        )

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

        if genome not in SUPPORTED_GENOMES:
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

        file_path = Path(file_record["path"])
        if not any(_is_within(file_path, d) for d in _managed_dirs(app)):
            return jsonify({"error": "Access denied"}), 403

        file_path = file_path.resolve()
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

        # Clean up files on disk. Only unlink paths inside the directories this
        # app owns, so a malformed or tampered DB row cannot delete elsewhere.
        managed = _managed_dirs(app)
        for f in db.get_all_files(analysis_id):
            file_path = Path(f["path"])
            if not any(_is_within(file_path, d) for d in managed):
                logger.warning(f"Refusing to delete file outside managed dirs: {file_path}")
                continue
            file_path.unlink(missing_ok=True)

        # Clean up result directory
        result_dir = Path(app.config["RESULTS_FOLDER"]) / str(analysis_id)
        if _is_within(result_dir, Path(app.config["RESULTS_FOLDER"])) and result_dir.exists():
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


def _recover_interrupted_runs(results_folder, db):
    """Reconcile output directories left behind by a process that died mid-swap.

    A run stages its output, swaps it into place, then commits. A crash in
    between can leave a `.<id>.incoming` (never swapped) or a `.<id>.previous`
    (swapped, commit unknown). publish_run sets status='complete' in the same
    transaction as the rows, so that flag says whether the commit landed:
    complete means the new output is live and the backup is stale, anything
    else means the swap must be undone to match the rows still in the database.

    Returns (restored, discarded) counts.
    """
    results_folder = Path(results_folder)
    restored = discarded = 0

    for entry in results_folder.glob(".*"):
        if not entry.is_dir():
            continue

        name = entry.name
        if name.endswith(".incoming"):
            # Never swapped in, so nothing references it.
            shutil.rmtree(entry, ignore_errors=True)
            discarded += 1
            continue

        if not name.endswith(".previous"):
            continue

        analysis_id = name[1 : -len(".previous")]
        analysis = db.get_analysis(int(analysis_id)) if analysis_id.isdigit() else None
        final_dir = results_folder / analysis_id

        if analysis is not None and analysis.get("status") == "complete":
            shutil.rmtree(entry, ignore_errors=True)
            discarded += 1
        else:
            shutil.rmtree(final_dir, ignore_errors=True)
            entry.rename(final_dir)
            restored += 1

    return restored, discarded


def _rebase(path, staging_dir, final_dir):
    """Re-point a path produced inside the staging directory at its final home.

    Preserves any subdirectory structure, so an output written to a nested path
    is registered where it will actually live rather than at the top level.
    """
    path = Path(path)
    try:
        relative = path.relative_to(staging_dir)
    except ValueError:
        # Not under staging, so the swap does not move it and its current
        # location stays valid. Rewriting it would point at a file that the
        # run never puts there.
        return str(path)
    return str(Path(final_dir) / relative)


def _swap_into_place(staging_dir, final_dir, analysis_id):
    """Move a completed run's output into final_dir, returning the kept-aside previous dir.

    Returns None when there was no previous output. The previous directory is
    retained so the caller can restore it if the database write fails. If the
    move itself fails the previous output is put back before raising, so the
    caller never sees a half-applied swap.
    """
    backup_dir = None
    if final_dir.exists():
        backup_dir = final_dir.parent / f".{analysis_id}.previous"
        shutil.rmtree(backup_dir, ignore_errors=True)
        final_dir.rename(backup_dir)
    try:
        staging_dir.rename(final_dir)
    except OSError:
        if backup_dir is not None:
            backup_dir.rename(final_dir)
        raise
    return backup_dir


def _restore_previous(backup_dir, final_dir):
    """Undo _swap_into_place after a failed publish."""
    if backup_dir is None:
        shutil.rmtree(final_dir, ignore_errors=True)
        return
    shutil.rmtree(final_dir, ignore_errors=True)
    backup_dir.rename(final_dir)


def run_analysis(analysis_id, db, socketio, results_folder=None):
    """Run analysis in background thread.

    Args:
        analysis_id: Analysis ID
        db: Database manager
        socketio: SocketIO instance for progress updates
        results_folder: Root directory for analysis output. Defaults to
            ~/.mutagene/results
    """
    from pathlib import Path

    from .analysis import run_cohort_analysis

    if results_folder is None:
        results_folder = Path.home() / ".mutagene" / "results"

    # Bound before the try so the finally can clean up even if setup fails.
    staging_dir = None
    # True once the run's results are committed; from then on a failure must not
    # relabel the analysis as failed.
    published = False

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
        final_dir = Path(results_folder) / str(analysis_id)
        # Write this run into a staging directory and move it into place only
        # once the analysis succeeds. Writing straight into final_dir let a
        # failed re-run overwrite the files that the still-published previous
        # results point at.
        staging_dir = Path(results_folder) / f".{analysis_id}.incoming"
        shutil.rmtree(staging_dir, ignore_errors=True)
        staging_dir.mkdir(parents=True, exist_ok=True)
        output_dir = staging_dir

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

        # Build the rows before touching the filesystem: malformed results must
        # fail while the previous run is still fully in place.
        file_rows = [
            (
                file_info["type"],
                Path(file_info["path"]).name,
                _rebase(file_info["path"], staging_dir, final_dir),
            )
            for file_info in results["files"]
        ]
        result_rows = [
            (result_type, results[result_type_key], None)
            for result_type, result_type_key in (
                ("profile", "profiles"),
                ("signatures", "signatures"),
                ("classification", "classification"),
                ("drivers", "drivers"),
                ("clustering", "clustering"),
                ("genome_warning", "genome_warning"),
            )
            if results.get(result_type_key)
        ]
        samples, mutations = results["samples"], results["mutations"]

        # The run succeeded, so move its output into place. The previous
        # directory is kept aside until the database write commits, so a failure
        # below can restore both the files and the rows together.
        backup_dir = _swap_into_place(staging_dir, final_dir, analysis_id)

        try:
            removed_results, removed_files = db.publish_run_results(
                analysis_id,
                samples,
                mutations,
                file_rows,
                result_rows,
            )
        except Exception:
            # Nothing was committed, so put the previous output back to match.
            _restore_previous(backup_dir, final_dir)
            raise

        # Past this point the run is committed and the analysis is 'complete'.
        # Anything that fails below is incidental and must not relabel it as
        # failed: doing so would contradict the committed rows and make startup
        # recovery restore the previous output over them.
        published = True

        if backup_dir is not None:
            shutil.rmtree(backup_dir, ignore_errors=True)
            if backup_dir.exists():
                logger.warning(
                    f"Could not remove the previous output directory {backup_dir}; "
                    f"it is no longer referenced and can be deleted"
                )

        if removed_results or removed_files:
            logger.info(
                f"Re-run of analysis {analysis_id} replaced {removed_results} previous "
                f"result(s) and {removed_files} output file record(s)"
            )

        socketio.emit("complete", {"analysis_id": analysis_id})

    except Exception as e:
        import traceback

        error_msg = f"{str(e)}\n{traceback.format_exc()}"
        if published:
            # The results are committed; reporting them is what broke. Telling
            # the client it failed would contradict what is stored.
            logger.exception(
                f"Analysis {analysis_id} completed and was stored, but reporting it failed"
            )
        else:
            db.update_analysis_status(analysis_id, "error", error_msg)
            socketio.emit("error", {"analysis_id": analysis_id, "error": str(e)})

    finally:
        # A failed run leaves its half-written output in staging; drop it rather
        # than let it accumulate or be picked up by the next attempt.
        if staging_dir is not None:
            shutil.rmtree(staging_dir, ignore_errors=True)


def start_server(host="127.0.0.1", port=5000, debug=False, open_browser=True):
    """Start the Flask development server.

    Args:
        host: Host to bind to
        port: Port to bind to
        debug: Enable debug mode
        open_browser: Automatically open browser

    Raises:
        ValueError: If debug mode is requested on a non-loopback bind address.
    """
    is_loopback = host in ("127.0.0.1", "localhost", "::1")

    if debug and not is_loopback:
        # Werkzeug's debug mode exposes an interactive console that executes
        # arbitrary Python for anyone who can reach it.
        raise ValueError(
            f"Refusing to enable debug mode while binding to {host}: the Werkzeug "
            "debugger allows remote code execution. Use --host 127.0.0.1 with --debug."
        )

    if not is_loopback:
        logger.warning(
            f"Binding to {host} exposes this server on the network. It has no "
            "authentication and runs analyses on uploaded files."
        )

    allowed_origins = [f"http://localhost:{port}", f"http://127.0.0.1:{port}"]
    allowed_hosts = ["localhost", "127.0.0.1", "::1"]
    if not is_loopback:
        # The hostname clients will use is unknown for a wildcard/LAN bind, so
        # the Host check cannot be applied; the warning above covers the risk.
        allowed_origins.append(f"http://{host}:{port}")
        allowed_hosts = None

    app, socketio = create_app(
        config={"ALLOWED_ORIGINS": allowed_origins, "ALLOWED_HOSTS": allowed_hosts}
    )

    print(f"""
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
""")

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
            # Flask-SocketIO refuses to start on Werkzeug unless this is set when
            # debug is off, which made `mutagene serve` exit immediately. This is
            # a single-user tool on loopback, which is what Werkzeug is fine for;
            # non-loopback binds are warned about above.
            allow_unsafe_werkzeug=True,
        )
    except KeyboardInterrupt:
        print("\n\nShutting down server...")
        sys.exit(0)
