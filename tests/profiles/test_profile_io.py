"""Reading profile files, and building profiles from counts.

The reader accepted anything that parsed as a float -- negative, NaN, inf -- and
let a repeated channel silently overwrite the earlier value. Both produce a
profile that looks ordinary and decomposes into nonsense.
"""

import math

import pytest

from mutagene.io.profile import get_profile_attributes_dict, read_profile_file
from mutagene.profiles.profile import get_mutational_profile


def write(tmp_path, text):
    path = tmp_path / "sample.profile"
    path.write_text(text)
    return str(path)


def full_profile(value="1.0", **overrides):
    lines = []
    for attr in get_profile_attributes_dict():
        p5, p3 = attr["context"]
        x, y = attr["mutation"]
        label = f"{p5}[{x}>{y}]{p3}"
        lines.append(f"{label}\t{overrides.get(label, value)}")
    return "\n".join(lines) + "\n"


class TestReadProfileFile:
    def test_failure_is_falsy(self):
        """(None, None) is truthy, so callers testing `if profile:` accepted it."""
        from mutagene.io.profile import read_profile_str

        assert not read_profile_str("nonsense\n")

    def test_a_well_formed_profile_reads(self, tmp_path):
        values = read_profile_file(write(tmp_path, full_profile()))

        assert values is not None
        assert len(values) == 96
        assert all(v == 1.0 for v in values)

    @pytest.mark.parametrize("bad", ["-1.0", "nan", "inf", "-inf"])
    def test_unusable_values_are_refused(self, tmp_path, bad):
        text = full_profile(**{"A[C>A]A": bad})

        assert read_profile_file(write(tmp_path, text)) is None

    def test_a_repeated_channel_is_refused(self, tmp_path):
        text = full_profile() + "A[C>A]A\t99.0\n"

        assert read_profile_file(write(tmp_path, text)) is None

    def test_zero_is_a_legitimate_count(self, tmp_path):
        values = read_profile_file(write(tmp_path, full_profile(value="0.0")))

        assert values is not None
        assert sum(values) == 0.0

    def test_a_missing_channel_defaults_to_zero(self, tmp_path):
        text = "\n".join(full_profile().splitlines()[1:]) + "\n"
        values = read_profile_file(write(tmp_path, text))

        assert values is not None
        assert len(values) == 96


class TestGetMutationalProfile:
    def test_counts_pass_through(self):
        values = get_mutational_profile({"AACT": 5.0}, counts=True)

        assert len(values) == 96
        assert sum(values) == 5.0

    def test_frequencies_sum_to_one(self):
        values = get_mutational_profile({"AACT": 3.0, "AACA": 1.0}, counts=False)

        assert sum(values) == pytest.approx(1.0)
        assert all(math.isfinite(v) for v in values)

    def test_an_empty_profile_says_so_instead_of_dividing_by_zero(self):
        with pytest.raises(ValueError, match="zero mutations"):
            get_mutational_profile({}, counts=False)

    def test_an_empty_profile_of_counts_is_all_zeros(self):
        values = get_mutational_profile({}, counts=True)

        assert values == [0] * 96
