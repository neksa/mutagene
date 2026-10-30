# Review the application-note draft

The LaTeX source and PDF describe local implementation and controlled benchmark results.
The draft does not claim a public release or clinical validation.
The draft covers the initial benchmark and does not yet include the sparse-selection follow-up.
The [sparse-selection report](../analysis/signature-benchmark-2026/sparse-report.html) documents that separate evaluation.

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
