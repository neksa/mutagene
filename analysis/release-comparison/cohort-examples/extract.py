"""Extract five candidate SBS96 signatures from the labelled PCAWG count matrix."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity

HERE = Path(__file__).resolve().parent


def extract(matrix, output):
    counts = pd.read_csv(matrix, sep="\t", index_col=0)
    catalog = pd.read_csv(HERE.parent / "performance/reference.tsv", sep="\t", index_col=0)
    if counts.shape != (96, 146) or not counts.index.equals(catalog.index):
        raise ValueError("Expected 146 PCAWG profiles in reference channel order")
    if (counts < 0).any().any() or not (counts.sum() > 0).all() or counts.isna().any().any():
        raise ValueError("Each tumor requires complete counts and positive total burden")
    profiles = counts.div(counts.sum()).T
    best_error = np.inf
    for seed in range(10):
        model = NMF(n_components=5, init="random", random_state=seed, tol=1e-6, max_iter=5000)
        exposures = model.fit_transform(profiles)
        if model.n_iter_ == model.max_iter:
            raise RuntimeError(f"NMF reaches the iteration limit for seed {seed}")
        if model.reconstruction_err_ < best_error:
            best_error = model.reconstruction_err_
            components = model.components_.copy()
            activities = exposures * components.sum(axis=1)

    order = np.argsort(-activities.mean(axis=0))
    components = components[order]
    components /= components.sum(axis=1, keepdims=True)
    names = [f"Component {i}" for i in range(1, len(order) + 1)]
    signatures = pd.DataFrame(components.T, index=counts.index, columns=names)
    similarity = cosine_similarity(signatures.T, catalog.T)
    matches = pd.DataFrame(
        {
            "component": names,
            "nearest_reference": catalog.columns[similarity.argmax(axis=1)],
            "cosine_similarity": similarity.max(axis=1),
        }
    )
    signatures.to_csv(output / "extracted-signatures.tsv", sep="\t")
    matches.to_csv(output / "signature-matches.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    extract(args.matrix, args.output)
