"""Exercise the installed web interface in Chromium with local test data."""

import shutil
import tempfile
from importlib.resources import files
from pathlib import Path
from threading import Thread
from unittest.mock import patch

from playwright.sync_api import expect, sync_playwright
from werkzeug.serving import make_server

from mutagene.webapp.database import DatabaseManager
from mutagene.webapp.server import create_app


def exercise(browser, origin, sample):
    page = browser.new_page(viewport={"width": 390, "height": 844})
    errors = []
    external = []
    page.on("pageerror", lambda error: errors.append(str(error)))

    def local_requests(route):
        if route.request.url.startswith(origin + "/"):
            route.continue_()
        else:
            external.append(route.request.url)
            route.abort()

    page.route("**/*", local_requests)
    page.goto(origin)
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    page.get_by_label("Mutation file").set_input_files(sample)
    page.get_by_label("Signature set").select_option("COSMICv2")
    page.get_by_role("button", name="Start analysis").click()
    expect(page.get_by_role("heading", name="Analysis Results")).to_be_visible(timeout=90000)
    page.get_by_role("button", name="Profiles", exact=True).click()
    graph = page.locator("#profileHistogram")
    expect(graph.locator(".main-svg").first).to_be_visible()
    channels = graph.evaluate(
        """graph => graph.data.flatMap(trace => trace.hovertemplate.map((label, i) =>
            [trace.name, label.split('<br>')[0], trace.y[i]]))"""
    )
    assert len(channels) == 96
    assert sum(count for _, _, count in channels) == 24
    for substitution, channel, count in channels:
        assert f"[{substitution}]" in channel
        row = page.locator("#profiles tbody tr").filter(has_text=channel)
        assert int(row.locator("td").nth(1).text_content()) == count
    page.get_by_role("button", name="Clustering", exact=True).click()
    expect(page.locator("#clustering tbody tr")).to_have_count(3)
    page.get_by_role("button", name="Driver ranking", exact=True).click()
    expect(page.locator("#hotspots")).to_contain_text("no driver ranking")
    page.set_viewport_size({"width": 1280, "height": 800})
    page.get_by_role("link", name="History", exact=True).click()
    expect(page.get_by_text("Complete", exact=True)).to_be_visible()
    page.get_by_role("link", name="New analysis").click()
    page.get_by_label("Mutation file").set_input_files(sample)
    page.route("**/api/analyze/*", lambda route: route.fulfill(status=503, json={"error": "Busy"}))
    page.get_by_role("button", name="Start analysis", exact=True).click()
    expect(page.get_by_role("alert")).to_contain_text("Busy")
    expect(page.get_by_role("button", name="Start analysis", exact=True)).to_be_enabled()
    assert not errors, errors
    assert not external, external
    page.close()


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    demo = files("mutagene").joinpath("data", "demo")
    genomes = root / "genomes"
    genomes.mkdir()
    shutil.copyfile(str(demo.joinpath("reference.2bit")), genomes / "hg19.2bit")
    with (
        patch("mutagene.webapp.genome_manager.DEFAULT_GENOMES_DIR", genomes),
        patch("mutagene.webapp.server.DatabaseManager", return_value=DatabaseManager(root / "db")),
    ):
        app, socketio = create_app(
            {
                "AUTO_DOWNLOAD_GENOMES": False,
                "UPLOAD_FOLDER": root / "uploads",
                "RESULTS_FOLDER": root / "results",
            }
        )
        server = make_server("127.0.0.1", 0, app, threaded=True)
        origin = f"http://127.0.0.1:{server.server_port}"
        app.config["ALLOWED_ORIGINS"].append(origin)
        socketio.server.eio.cors_allowed_origins = [origin]
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch()
                exercise(browser, origin, str(demo.joinpath("sample.maf")))
                browser.close()
        finally:
            server.shutdown()
            thread.join()
            server.server_close()
            app.config["ANALYSIS_EXECUTOR"].shutdown()
print("Browser checks pass for upload, results, channel counts, history, and start failure.")
