"""Read and write labelled mutation matrices for Python and R workflows."""

from itertools import product

import numpy as np
import pandas as pd

from mutagene.io.profile import get_profile_attributes_dict


def channel_labels(context="SBS96"):
    """Return labels in MutaGene signature order."""
    if context == "SBS1536":
        flanks = ["".join(pair) for pair in product("ACGT", repeat=2)]
        return [
            f"{left}[{substitution}]{right}"
            for substitution in ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
            for left in flanks
            for right in flanks
        ]
    if context == "ID83":
        labels = [
            f"1:{event}:{base}:{repeat}"
            for event in ("Del", "Ins")
            for base in "CT"
            for repeat in range(6)
        ]
        labels.extend(
            f"{length}:{event}:R:{repeat}"
            for event in ("Del", "Ins")
            for length in range(2, 6)
            for repeat in range(6)
        )
        for length in range(2, 6):
            for homology in range(1, 6 if length == 5 else length):
                labels.append(f"{length}:Del:M:{homology}")
        return labels
    if context != "SBS96":
        raise ValueError(f"Unknown mutation context: {context}")
    return [
        f"{a['context'][0]}[{a['mutation'][0]}>{a['mutation'][1]}]{a['context'][1]}"
        for a in get_profile_attributes_dict()
    ]


def _read_matrix(stream, context):
    header = stream.readline().rstrip("\r\n").split("\t")
    if len(header) < 2 or header[0] not in {"Type", "MutationType"}:
        raise ValueError(
            "A matrix requires a Type or MutationType column and at least one data column"
        )
    if len(set(header)) != len(header) or any(not name.strip() for name in header):
        raise ValueError("Matrix column names must be nonempty and unique")
    frame = pd.read_csv(stream, sep="\t", header=None, names=header, index_col=0)
    labels = channel_labels(context)
    if not frame.index.is_unique or set(frame.index) != set(labels):
        raise ValueError(
            f"A matrix requires each of the {len(labels)} {context} channels exactly once"
        )
    return frame.loc[labels].astype(float)


def read_reference_matrix(stream, context):
    """Read a labelled channel-by-signature matrix and preserve signature names."""
    frame = _read_matrix(stream, context)
    values = frame.to_numpy()
    if not np.isfinite(values).all() or (values < 0).any() or (values.sum(axis=0) <= 0).any():
        raise ValueError("Reference columns require finite nonnegative values and positive mass")
    return values, list(frame.columns)


def read_count_matrix(stream, context="SBS96"):
    """Read a Type column and one count column per sample."""
    frame = _read_matrix(stream, context)
    values = frame.to_numpy()
    if not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(values != np.floor(values)):
        raise ValueError("Matrix values must be finite, nonnegative integer counts")
    if np.any(values.sum(axis=0) == 0):
        raise ValueError("Every matrix sample must contain at least one mutation")
    return {name: values[:, i] for i, name in enumerate(frame.columns)}


def write_count_matrix(stream, profiles, context="SBS96"):
    """Write one sample per column and preserve the canonical channel order."""
    pd.DataFrame(profiles, index=channel_labels(context)).rename_axis("Type").to_csv(
        stream, sep="\t", float_format="%.0f"
    )
