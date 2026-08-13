"""Resolution of the -g argument against ~/.mutagene/genomes (issue #93)."""

import subprocess
import sys

import pytest

from mutagene.cli.genome_resolver import (
    DEFAULT_GENOMES_DIR,
    GenomeNotFoundError,
    resolve_genome,
)


@pytest.fixture
def genomes_dir(tmp_path):
    """A stand-in for ~/.mutagene/genomes holding a downloaded hg19."""
    d = tmp_path / "genomes"
    d.mkdir()
    (d / "hg19.2bit").write_bytes(b"")
    return d


def test_bare_name_resolves_to_genomes_dir(genomes_dir):
    assert resolve_genome("hg19", genomes_dir) == str(genomes_dir / "hg19.2bit")


def test_bare_name_with_extension_resolves_to_genomes_dir(genomes_dir):
    assert resolve_genome("hg19.2bit", genomes_dir) == str(genomes_dir / "hg19.2bit")


def test_full_path_is_used_as_is(tmp_path, genomes_dir):
    custom = tmp_path / "custom" / "hg19.2bit"
    custom.parent.mkdir()
    custom.write_bytes(b"")
    assert resolve_genome(str(custom), genomes_dir) == str(custom)


def test_relative_path_is_used_as_is(tmp_path, genomes_dir, monkeypatch):
    (tmp_path / "local.2bit").write_bytes(b"")
    monkeypatch.chdir(tmp_path)
    assert resolve_genome("./local.2bit", genomes_dir) == "./local.2bit"


def test_current_directory_fallback(tmp_path, genomes_dir, monkeypatch):
    """A bare name still works from a directory holding <name>.2bit."""
    (tmp_path / "mm10.2bit").write_bytes(b"")
    monkeypatch.chdir(tmp_path)
    assert resolve_genome("mm10", genomes_dir) == "mm10"


def test_unknown_genome_names_both_locations(tmp_path, genomes_dir, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(GenomeNotFoundError) as excinfo:
        resolve_genome("nosuch", genomes_dir)

    message = str(excinfo.value)
    assert str(genomes_dir / "nosuch.2bit") in message
    assert str(tmp_path / "nosuch.2bit") in message
    assert "mutagene fetch genome -g nosuch" in message


def test_genomes_dir_defaults_to_mutagene_home(tmp_path, monkeypatch):
    """Without an explicit directory, ~/.mutagene/genomes is the one searched."""
    monkeypatch.chdir(tmp_path)
    with pytest.raises(GenomeNotFoundError) as excinfo:
        resolve_genome("not_a_genome")
    assert str(DEFAULT_GENOMES_DIR / "not_a_genome.2bit") in str(excinfo.value)


def test_default_genomes_dir_is_the_fetch_destination():
    assert DEFAULT_GENOMES_DIR.name == "genomes"
    assert DEFAULT_GENOMES_DIR.parent.name == ".mutagene"


def test_cli_imports_do_not_require_flask():
    """The [web] extra is optional: importing the CLI must not pull in Flask."""
    # Blocking the modules makes any import of them from the CLI path fail
    code = (
        "import sys;"
        "sys.modules['flask'] = None;"
        "sys.modules['flask_socketio'] = None;"
        "import mutagene.__main__;"
        "from mutagene.cli.genome_resolver import resolve_genome"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
