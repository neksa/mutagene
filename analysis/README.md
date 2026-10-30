# Review the research evidence

This directory contains historical benchmark sources, summary tables, reports, and provenance records.
The runtime PRs remain separate from this research draft.

## Choose the relevant assessment

- [Likelihood repair](signature-benchmark-2026/repaired-report.qmd) supersedes conclusions that depend on the original solver.
- [Background selection](signature-benchmark-2026/background-report.qmd) compares fixed, joint, and block selection.
- [Low-burden search](low-burden-2026/report.qmd) separates tuning, validation, and sensitivity panels.
- [Threshold scenarios](search-scenarios-2026/report.qmd) measures the discontinuity at 1,000 mutations.
- [signeR comparison](signature-benchmark-2026/signer-report.qmd) separates sampling completion from convergence evidence.
- [Original assessment](signature-benchmark-2026/report.qmd) preserves the initial analysis and its reproduction commands.

These synthetic experiments do not establish clinical validity or universal attribution superiority.
Historical reports use the default and source snapshots from their recorded experiment.
They do not describe the default on every PR branch.
The separate default-change PR controls adoption.

## Verify this source-and-summary snapshot

From the repository root, run:

```sh
python analysis/verify_publication.py
```

The publication manifest covers the files in this review snapshot.
It does not certify the underlying scientific claims.

## Reproduce or audit complete experiments

The original experiment manifests also name bulk data that this source-and-summary snapshot does not include.
`local-artifacts.json` records the omitted paths, byte counts, and SHA256 hashes.
These include per-fit output, simulation arrays, and third-party reference matrices.
The original local checkout retains those files.
No public data archive or configured DVC remote accompanies this PR.

Full arithmetic verifiers require the complete matching evidence bundle.
Some historical specifications contain absolute paths from the original workstation.
Do not interpret missing files or path errors as a successful verification.
Use the report's pinned source revision and environment to regenerate an experiment in a fresh output directory.
Do not overwrite historical results when testing current code.

Before public redistribution, review third-party data terms and select durable artifact storage.
Preserve source snapshots and original manifests when archiving results.
Current source edits can invalidate historical checks against live package paths.

## Review the manuscript separately

The [application-note draft](../manuscript/README.md) predates the numerical repair and sparse-selection assessment.
Its PDF does not incorporate those results.
Revise its claims and figures before submission.
The [extended-context proposal](../docs/extended-context-logos.md) describes a future feature.
