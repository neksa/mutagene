"""Learn one extra mutational signature while holding a supplied reference catalog fixed."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls


def fit(basis, profiles):
    q, r = np.linalg.qr(basis, mode="reduced")
    return np.column_stack([nnls(r, target)[0] for target in (q.T @ profiles).T])


def unit_columns(x):
    totals = x.sum(axis=0)
    if np.any(totals <= 0):
        raise ValueError("Each spectrum requires positive mass")
    return x / totals


def predict(basis, profiles):
    weights = fit(basis, profiles)
    # Equal reference weights give a prediction when no reference overlaps a fitting observation.
    weights[:, weights.sum(axis=0) == 0] = 1
    weights = unit_columns(weights)
    return basis @ weights, weights


def cosine(a, b):
    return (a / np.linalg.norm(a, axis=0)).T @ (b / np.linalg.norm(b, axis=0))


def learn(basis, training):
    baseline = basis @ fit(basis, training)
    residual = np.maximum(training - baseline, 0)
    if np.linalg.norm(residual) <= 1e-12 * np.linalg.norm(training):
        raise ValueError(
            "The supplied references explain the training profiles within numerical precision"
        )
    candidate = unit_columns(np.abs(np.linalg.svd(residual, full_matrices=False)[0][:, 0]))
    previous_loss = np.inf
    for _ in range(2000):
        joint = np.column_stack([basis, candidate])
        weights = fit(joint, training)
        loss = np.square(training - joint @ weights).sum()
        if abs(previous_loss - loss) <= 1e-6 * loss:
            return candidate
        extra = weights[-1]
        updated = unit_columns(np.maximum((training - basis @ weights[:-1]) @ extra, 0))
        previous_loss = loss
        candidate = updated
    raise RuntimeError("Candidate learning did not converge")


def compare_candidate(counts_path, reference_path, out, withhold=None):
    out = Path(out)
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    catalog = pd.read_csv(reference_path, sep="\t", index_col=0)
    if (
        not counts.index.is_unique
        or not catalog.index.is_unique
        or set(catalog.index) != set(counts.index)
    ):
        raise ValueError("Reference and count matrices must share unique channel labels")
    catalog = catalog.loc[counts.index]
    for frame in (catalog, counts):
        if frame.empty:
            raise ValueError("Supply at least one reference and at least one channel")
        if (
            not np.isfinite(frame).all().all()
            or (frame < 0).any().any()
            or (frame.sum() <= 0).any()
        ):
            raise ValueError("Each profile requires finite nonnegative values and positive mass")
    if (counts != np.floor(counts)).any().any() or counts.shape[1] < 6:
        raise ValueError("Supply integer counts for at least six tumors")
    if withhold and len(catalog.columns) < 2:
        raise ValueError("The control requires at least two reference signatures")
    if withhold and withhold not in catalog:
        raise ValueError("The withheld reference must appear in the catalog")
    rng = np.random.default_rng(317)
    histology = np.array([name.split("::")[0] if "::" in name else "cohort" for name in counts])
    test = np.zeros(counts.shape[1], dtype=bool)
    for group in sorted(set(histology)):
        members = np.flatnonzero(histology == group)
        test[rng.permutation(members)[: len(members) // 3]] = True
    if not test.any():
        raise ValueError("At least one group requires three tumors for a held-out comparison")
    train = unit_columns(counts.loc[:, ~test].to_numpy())
    held = counts.loc[:, test].to_numpy(dtype=np.int64)
    totals = held.sum(axis=0)
    if (totals < 2).any():
        raise ValueError("Each held-out tumor needs at least two mutations")
    fitted_half = rng.binomial(held, 0.5)
    empty = (fitted_half.sum(axis=0) == 0) | (fitted_half.sum(axis=0) == totals)
    while empty.any():
        fitted_half[:, empty] = rng.binomial(held[:, empty], 0.5)
        empty = (fitted_half.sum(axis=0) == 0) | (fitted_half.sum(axis=0) == totals)
    observations = unit_columns(held - fitted_half)
    fitted_half = unit_columns(fitted_half)
    rows, matches, channels, scores = [], [], {}, []
    experiments = [("Full catalog", list(catalog))]
    if withhold:
        experiments.append((f"Withhold {withhold}", list(catalog.drop(columns=withhold))))
    for experiment, names in experiments:
        basis = unit_columns(catalog[names].to_numpy())
        candidate = learn(basis, train)
        channels[experiment] = candidate
        refit = fit(basis, candidate[:, None])[:, 0]
        mixture = basis @ refit
        mixture_cosine = 0.0
        if np.linalg.norm(mixture) > 1e-12 * np.linalg.norm(candidate):
            refit = unit_columns(refit)
            mixture = unit_columns(mixture)
            mixture_cosine = cosine(candidate, mixture)
        else:
            refit[:] = np.nan
            mixture[:] = np.nan
        channels[experiment + " mixture"] = mixture
        similarities = cosine(candidate, basis)
        for name, similarity, weight in zip(names, similarities, refit):
            matches.append((experiment, name, similarity, weight))
        baseline, _ = predict(basis, fitted_half)
        augmented = np.column_stack([basis, candidate])
        prediction, weights = predict(augmented, fitted_half)
        before = np.abs(baseline - observations).sum(axis=0) / 2
        after = np.abs(prediction - observations).sum(axis=0) / 2
        for sample, a, b, w in zip(counts.columns[test], before, after, unit_columns(weights)[-1]):
            scores.append((experiment, sample, a, b, w))
        gain = before - after
        boot = np.random.default_rng(318).choice(gain, size=(5000, len(gain))).mean(axis=1)
        rows.append(
            (
                experiment,
                len(names),
                train.shape[1],
                len(gain),
                names[similarities.argmax()],
                similarities.max(),
                mixture_cosine,
                (
                    cosine(candidate, catalog[withhold].values)
                    if withhold and withhold not in names
                    else np.nan
                ),
                gain.mean(),
                *np.quantile(boot, [0.025, 0.975]),
                (gain > 1e-12).sum(),
            )
        )
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(channels, index=counts.index).to_csv(out / "candidate-profiles.tsv", sep="\t")
    pd.DataFrame(matches, columns=["experiment", "reference", "cosine", "mixture_weight"]).to_csv(
        out / "candidate-matches.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        scores,
        columns=["experiment", "tumor", "reference_error", "augmented_error", "candidate_exposure"],
    ).to_csv(out / "held-out.tsv", sep="\t", index=False)
    summary = pd.DataFrame(
        rows,
        columns=[
            "experiment",
            "references",
            "training_tumors",
            "test_tumors",
            "nearest",
            "cosine",
            "mixture_cosine",
            "withheld_reference_cosine",
            "mean_gain",
            "ci_low",
            "ci_high",
            "improved_tumors",
        ],
    )
    summary.to_csv(out / "candidate-summary.tsv", index=False, sep="\t")
    return summary
