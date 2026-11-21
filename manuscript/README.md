# Review the application-note draft

The [current preprint](preprint.pdf) covers attribution benchmarks, motif enrichment, candidate signature extraction, PCAWG exposures, and published treatment outcomes.
The [extended-context examples](../analysis/release-comparison/extended-contexts/README.md) demonstrate SBS1536 and ID83 fitting.
The [matched PCAWG analysis](../analysis/release-comparison/pcawg-penta/README.md) compares SBS1536 candidates with supplied signatures.
It also reports motif enrichment and BGI driver rankings.
The [external comparisons](../analysis/release-comparison/cohort-examples/replication.tsv) assess all five PCAWG components.
See the [benchmark tables](../analysis/release-comparison/performance/README.md) and [cohort figure sources](../analysis/release-comparison/cohort-examples/README.md).
To rebuild it, run `pdflatex -halt-on-error preprint.tex` twice from this directory.

## Historical draft

This historical draft predates the likelihood repair and sparse-selection assessment.
Do not use its claims to describe the repaired implementation.
The LaTeX source and PDF preserve the initial controlled benchmark results.
The draft does not claim a public release or clinical validation.
The draft covers the initial benchmark and does not yet include the sparse-selection follow-up.
The [sparse-selection report](../analysis/signature-benchmark-2026/sparse-report.html) documents that separate evaluation.
The [likelihood repair report](../analysis/signature-benchmark-2026/repaired-report.html) documents later numerical corrections and reassessment.
The historical PDF preserves those earlier results.

## Build

Run these commands from this directory with a LaTeX distribution on PATH:

```sh
pdflatex -interaction=nonstopmode -halt-on-error application-note.tex
pdflatex -interaction=nonstopmode -halt-on-error application-note.tex
```

The figure comes from `../analysis/signature-benchmark-2026/figures/attribution.pdf`.
The Quarto report provides detailed methods and the evidence bundle.

## Before submission

Confirm authorship, affiliations, funding, competing interests, and acknowledgments.
Review the scientific claims and select an archival software release.
Confirm third-party data redistribution terms and final data-availability wording.
The working draft contains no invented author list or funding statement.
