"""Generate fresh empirical-mixture panels across low, intermediate, and high burdens."""

from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import CANCERS, CONDITIONS, OUT_SIGNATURES, channel_order

from mutagene.signatures import get_dummy_signatures_lists


def generate_panel(seed, samples, burdens, cutoff=True):
    source = Path("/tmp/mutagene-sigfittest-2024")
    reference = pd.read_csv(source / "input/COSMIC_v3_SBS_GRCh38.txt", sep="\t", index_col=0)
    generating = pd.read_csv(source / "input/COSMIC_v3.3.1_SBS_GRCh38.txt", sep="\t", index_col=0)
    reference, generating = (x.loc[channel_order()] for x in (reference, generating))
    reference /= reference.sum()
    generating /= generating.sum()
    dummy = np.array([v for _, v in get_dummy_signatures_lists()]).T
    rows, catalogs, truths = [], [], []
    for ci, cancer in enumerate(CANCERS):
        empirical = pd.read_csv(
            source / f"cosmic tissue data/signature_contributions_WGS_{cancer}.dat",
            sep="\t",
            index_col=0,
        )
        rng = np.random.default_rng(np.random.SeedSequence([seed, ci]))
        selected = rng.choice(len(empirical), samples, replace=True)
        for sample_index, patient in enumerate(selected):
            original = empirical.iloc[patient].to_numpy(float)
            original /= original.sum()
            for burden in burdens:
                weights = original.copy()
                if cutoff:
                    weights[weights * burden < 10] = 0
                if weights.sum() == 0:
                    raise ValueError("The truth cutoff removes every component")
                weights /= weights.sum()
                for condition in CONDITIONS:
                    truth = pd.Series(0.0, index=generating.columns)
                    truth.loc[empirical.columns] = weights
                    null_truth = 0.0 if condition == "clean" else 0.2
                    truth *= 1 - null_truth
                    if condition == "out_of_reference_20":
                        truth.loc[rng.choice(OUT_SIGNATURES, 2, replace=False)] = 0.1
                    probabilities = generating.to_numpy() @ truth.to_numpy()
                    if condition == "uniform_20":
                        probabilities += 0.2 / 96
                    elif condition == "substitution_noise_20":
                        probabilities += 0.2 * dummy @ np.array([0.05, 0.1, 0.6, 0.1, 0.1, 0.05])
                    probabilities /= probabilities.sum()
                    counts = rng.multinomial(burden, probabilities)
                    rows.append(
                        dict(
                            sample=f"{seed}_{cancer}_{sample_index}_{burden}_{condition}",
                            cancer=cancer,
                            sample_index=sample_index,
                            empirical_index=int(patient),
                            burden=burden,
                            condition=condition,
                            null_truth=null_truth,
                        )
                    )
                    catalogs.append(counts)
                    truths.append(truth.to_numpy())
    return dict(
        counts=np.array(catalogs),
        truth=np.array(truths),
        reference=reference.to_numpy(),
        names=reference.columns.to_numpy(dtype=str),
        generating_names=generating.columns.to_numpy(dtype=str),
    ), pd.DataFrame(rows)
