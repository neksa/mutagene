"""Matrix interchange must preserve channel labels and count semantics."""

import io

import numpy as np
import pytest

from mutagene.io.matrix import channel_labels, read_count_matrix, write_count_matrix


def test_reorders_channels_by_label():
    text = "Type\tA\n" + "".join(
        f"{c}\t{i + 1}\n" for i, c in reversed(list(enumerate(channel_labels())))
    )
    profiles = read_count_matrix(io.StringIO(text))
    np.testing.assert_array_equal(profiles["A"], np.arange(1, 97))


@pytest.mark.parametrize("bad_value", ["NaN", "inf", "-1", "0.5"])
def test_rejects_invalid_counts(bad_value):
    text = "Type\tA\n" + "".join(
        f"{c}\t{bad_value if i == 0 else 1}\n" for i, c in enumerate(channel_labels())
    )
    with pytest.raises(ValueError, match="integer counts"):
        read_count_matrix(io.StringIO(text))


def test_roundtrip_preserves_each_sample():
    profiles = {"A": np.arange(96), "B": np.arange(96)[::-1]}
    stream = io.StringIO()
    write_count_matrix(stream, profiles)
    stream.seek(0)
    result = read_count_matrix(stream)
    for name in profiles:
        np.testing.assert_array_equal(result[name], profiles[name])


def test_rejects_duplicate_sample_names():
    with pytest.raises(ValueError, match="unique"):
        read_count_matrix(io.StringIO("Type\tA\tA\n"))
