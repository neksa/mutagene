"""The benchmark subcommand (GitHub issue #114).

Every mode failed. Pairwise modes cast the signature set to an int before a
lookup keyed on strings; multiple modes used variables only ever assigned in the
pairwise branch; and the module called three library functions whose signatures
had changed underneath it.
"""

import glob
import os

import pytest

from mutagene.benchmark.generate_benchmark import (
    aggregate_benchmarks,
    gen_sample_2combinations,
    run_benchmark_2combinations,
)
from mutagene.benchmark.multiple_benchmark import multiple_benchmark, multiple_benchmark_run
from mutagene.io.profile import read_signatures

# MGA is the smallest set, so the tests stay quick.
SIGNATURE_SET = "5"


@pytest.fixture(scope="module")
def signatures():
    return read_signatures(SIGNATURE_SET)


class TestSignatureSetLookup:
    def test_the_set_name_is_a_string(self):
        """int("5") is not a key; this is what broke every pairwise mode."""
        W, names = read_signatures(SIGNATURE_SET)
        assert W.shape[0] == 96 and len(names) == 5

    def test_an_int_is_rejected(self):
        with pytest.raises(ValueError, match="Unknown signature set"):
            read_signatures(5)


class TestPairwise:
    def test_generate_then_run_then_aggregate(self, tmp_path, signatures):
        W, names = signatures
        root = str(tmp_path)

        gen_sample_2combinations(root, names, W, ratio=0.7, noise_level=0.1, n_mutations=50)
        profiles = glob.glob(os.path.join(root, "**", "*.profile"), recursive=True)
        assert profiles, "generation produced nothing"

        run_benchmark_2combinations(root, SIGNATURE_SET, names, W, force=True)
        infos = glob.glob(os.path.join(root, "**", "*.info"), recursive=True)
        assert infos, "the run produced no decompositions"

        aggregate_benchmarks(root)


class TestMultiple:
    def test_generate_then_run(self, tmp_path):
        root = str(tmp_path)

        multiple_benchmark(
            data_root=root, signature_sets=[SIGNATURE_SET], replicates=2, processes=1
        )
        profiles = glob.glob(os.path.join(root, "multiple", "*.profile"))
        assert len(profiles) == 2, f"expected 2 generated profiles, got {len(profiles)}"

        W, names = read_signatures(SIGNATURE_SET)
        multiple_benchmark_run(SIGNATURE_SET, names, W, force=True, data_root=root, processes=1)
        assert glob.glob(os.path.join(root, "multiple", "*.info"))

    def test_run_without_generated_data_warns_rather_than_failing(self, tmp_path, caplog):
        import logging

        W, names = read_signatures(SIGNATURE_SET)
        os.makedirs(os.path.join(str(tmp_path), "multiple"), exist_ok=True)

        with caplog.at_level(logging.WARNING):
            multiple_benchmark_run(SIGNATURE_SET, names, W, data_root=str(tmp_path), processes=1)

        assert any("No generated profiles" in r.message for r in caplog.records)

    def test_the_root_argument_is_honoured(self, tmp_path):
        """The multiple modes wrote to a hardcoded data/benchmark/multiple."""
        root = str(tmp_path / "custom")
        before = sorted(glob.glob("data/benchmark/multiple/*"))

        multiple_benchmark(
            data_root=root, signature_sets=[SIGNATURE_SET], replicates=1, processes=1
        )

        assert glob.glob(os.path.join(root, "multiple", "*.profile"))
        assert before == sorted(
            glob.glob("data/benchmark/multiple/*")
        ), "the hardcoded default directory was written to"
