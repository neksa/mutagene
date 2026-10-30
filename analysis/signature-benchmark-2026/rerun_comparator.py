"""Repeat SPA with its required channel order and preserve frozen MutaGene fits."""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from run_benchmark import score, sha256
from SigProfilerAssignment import Analyzer
from sigProfilerPlotting import get_context_reference

root = Path(__file__).resolve().parent / "results"
data = np.load(root / "inputs.npz")
metadata = pd.read_csv(root / "samples.tsv", sep="\t")
reference = pd.read_csv(root / "reference.tsv", sep="\t", index_col=0)
reference.loc[get_context_reference("96")].to_csv(
    root / "spa_reference.tsv", sep="\t", index_label="MutationType"
)
start = time.perf_counter()
Analyzer.cosmic_fit(
    samples=str(root / "catalogs.tsv"),
    output=str(root / "spa_corrected"),
    signature_database=str(root / "spa_reference.tsv"),
    input_type="matrix",
    genome_build="GRCh38",
    cosmic_version=3,
    make_plots=False,
    sample_reconstruction_plots=False,
    export_probabilities=False,
    cpu=4,
)
elapsed = time.perf_counter() - start
activities = pd.read_csv(
    root / "spa_corrected/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    sep="\t",
    index_col=0,
)
metrics = pd.read_csv(root / "metrics.tsv", sep="\t")
keep = metrics.method != "sigprofilerassignment"
estimates = np.load(root / "estimates.npz")
knowns = list(estimates["known"][keep])
nulls = list(estimates["null"][keep])
rows = []
for i, row in metadata.iterrows():
    known = activities.loc[row["sample"]].reindex(reference.columns).to_numpy(float) / row["burden"]
    null = np.zeros(6)
    result = row.to_dict()
    result.update(method="sigprofilerassignment", seconds=elapsed / len(metadata))
    result.update(
        score(
            known,
            null,
            data["truth"][i],
            data["counts"][i],
            reference.columns.tolist(),
            data["generating_names"].tolist(),
            data["reference"],
            row["null_truth"],
        )
    )
    rows.append(result)
    knowns.append(known)
    nulls.append(null)
metrics = pd.concat([metrics[keep], pd.DataFrame(rows)], ignore_index=True)
metrics.to_csv(root / "metrics.tsv", sep="\t", index=False)
np.savez_compressed(root / "estimates.npz", known=np.array(knowns), null=np.array(nulls))
metrics.groupby(["condition", "burden", "method"]).mean(numeric_only=True).to_csv(
    root / "summary.tsv", sep="\t"
)
provenance = json.loads((root / "provenance.json").read_text())
provenance["initial_runner_sha256"] = provenance["runner_sha256"]
provenance["runner_sha256"] = sha256(Path(__file__).with_name("run_benchmark.py"))
provenance["comparator_rerun_sha256"] = sha256(__file__)
provenance["comparator_elapsed_seconds"] = elapsed
provenance["notes"].append(
    "SPA repeats with its canonical reference channel order; frozen MutaGene fits remain identical."
)
provenance["output_hashes"] = {
    p.name: sha256(p) for p in root.iterdir() if p.is_file() and p.name != "provenance.json"
}
(root / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print("Corrected comparator complete", flush=True)
