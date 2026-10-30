"""Read and write labelled SBS96 count matrices for Python and R workflows."""

import numpy as np
import pandas as pd

from mutagene.io.profile import get_profile_attributes_dict


def channel_labels():
    """Return labels in MutaGene signature order."""
    return [
        f"{a['context'][0]}[{a['mutation'][0]}>{a['mutation'][1]}]{a['context'][1]}"
        for a in get_profile_attributes_dict()
    ]


def read_count_matrix(stream):
    """Read a Type column and one count column per sample."""
    # Inspect the header before pandas renames duplicate column names.
    header = stream.readline().rstrip("\r\n").split("\t")
    if len(header) < 2 or header[0] != "Type":
        raise ValueError("A count matrix requires a Type column and at least one sample column")
    if len(set(header[1:])) != len(header[1:]) or any(not name for name in header[1:]):
        raise ValueError("Matrix sample names must be nonempty and unique")
    frame = pd.read_csv(stream, sep="\t", header=None, names=header, index_col=0)
    labels = channel_labels()
    if not frame.index.is_unique or set(frame.index) != set(labels):
        raise ValueError("A count matrix requires each of the 96 SBS channels exactly once")
    values = frame.loc[labels].to_numpy(dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(values != np.floor(values)):
        raise ValueError("Matrix values must be finite, nonnegative integer counts")
    if np.any(values.sum(axis=0) == 0):
        raise ValueError("Every matrix sample must contain at least one mutation")
    return {name: values[:, i] for i, name in enumerate(frame.columns)}


def write_count_matrix(stream, profiles):
    """Write one sample per column and preserve the canonical channel order."""
    pd.DataFrame(profiles, index=channel_labels()).rename_axis("Type").to_csv(
        stream, sep="\t", float_format="%.0f"
    )
