"""mutagene.reports.nci60 (GitHub issue #115).

The module star-imported three modules that are not in the package, inside a
`try: ... except ImportError: pass`. It therefore imported cleanly and raised
NameError on seven names as soon as it did any work. Nothing caught it: star
imports stop ruff tracking undefined names, and nothing imported the module.
"""

import ast
import builtins
import inspect
import os

import mutagene.reports.nci60 as nci60


def test_the_module_imports():
    assert inspect.ismodule(nci60)


def test_no_name_is_used_without_being_defined():
    """The check that would have caught the original breakage."""
    source = inspect.getsource(nci60)
    tree = ast.parse(source)

    defined = set(dir(builtins))
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            defined.add(node.name)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            defined.add(node.id)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                if alias.name != "*":
                    defined.add(alias.asname or alias.name.split(".")[0])
        elif isinstance(node, ast.arg):
            defined.add(node.arg)
        elif isinstance(node, ast.ExceptHandler) and node.name:
            defined.add(node.name)

    used = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}

    assert not (used - defined), f"undefined at runtime: {sorted(used - defined)}"


def test_no_star_imports_remain():
    """They are what let the undefined names go unnoticed."""
    tree = ast.parse(inspect.getsource(nci60))
    starred = [
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
        for alias in node.names
        if alias.name == "*"
    ]
    assert not starred, f"star imports: {starred}"


def test_the_expected_entry_points_exist():
    for name in ("make_NCI60_report", "analyze_nci60_samples", "export_maf"):
        assert callable(getattr(nci60, name))


class TestPaths:
    def test_the_genome_path_is_not_machine_specific(self):
        assert not nci60.TWOBIT_GENOMES_PATH.startswith("/net/")
        assert "gonceare" not in nci60.TWOBIT_GENOMES_PATH

    def test_the_genome_path_is_overridable(self, monkeypatch):
        import importlib

        monkeypatch.setenv("MUTAGENE_GENOMES", "/tmp/genomes")
        importlib.reload(nci60)
        try:
            assert nci60.TWOBIT_GENOMES_PATH == "/tmp/genomes"
        finally:
            monkeypatch.delenv("MUTAGENE_GENOMES")
            importlib.reload(nci60)

    def test_no_home_directory_is_hardcoded_anywhere(self):
        source = inspect.getsource(nci60)
        for marker in ("/Users/", "/net/pan1", "anaconda3"):
            assert marker not in source, f"{marker} is still hardcoded"


def test_reports_are_reachable_from_the_package():
    """Nothing imported this module, which is why the breakage stayed hidden."""
    assert os.path.exists(os.path.join("mutagene", "reports", "__init__.py"))
