# Read or rebuild the application note

The [preprint](preprint.pdf) follows the analysis from tumor profiles to driver candidates.
Pages 2–4 cover SBS96, SBS1536, and ID83 fitting.
Page 5 shows motif profiles, enrichment, and bladder tumor stage.
Page 6 compares a learned signature with supplied references.
Pages 7–8 compare hotspot recovery against recurrence and examine gene-level ranking.

## Reproduce all seven figures and the PDF

Use Python 3.12, DVC 3.66.1, and a LaTeX installation with `pdflatex` on PATH.
The LaTeX source uses `booktabs`, `geometry`, `hyperref`, `lmodern`, and `microtype`.
From the repository root, run:

```sh
python3.12 -m venv .venv-paper
. .venv-paper/bin/activate
python -m pip install -e . -r manuscript/requirements.txt 'dvc==3.66.1'
dvc repro --force paper
dvc status dvc.yaml
```

The [DVC pipeline](../dvc.yaml) generates seven vector figures and the hotspot recovery table before two LaTeX passes.
The [lockfile](../dvc.lock) records checksums. Git stores the compact inputs and outputs. The build needs no DVC remote.
A fixed rendering date keeps PDF metadata stable. Font and TeX versions can affect binary output across systems.

## Evidence behind the figures

| Figure | Starting data and computation |
|---|---|
| Attribution accuracy | [Matched tool scores](../analysis/release-comparison/performance/README.md) from the archived R and Python benchmarks. |
| SBS96 exposures and outcomes | [PCAWG weights, TOPACIO event counts, and PETREMAC response counts](../analysis/release-comparison/cohort-examples/README.md). |
| SBS1536 exposures | [Matched PCAWG counts, fitted exposures, and fit checks](../analysis/release-comparison/pcawg-penta/README.md). |
| Motifs | The pipeline computes seven motif fractions per tumor, 7,715 histology comparisons, simulation controls, and the bladder stage comparison. |
| Candidate fingerprints | [Held-out scores and learned native-resolution profiles](../analysis/release-comparison/pcawg-penta/README.md). The plot aligns all channels by label. |
| Driver recovery | The pipeline fits the pooled bladder spectrum and scores the same missense changes under empirical and reconstructed backgrounds. |
| Extended contexts | [Published SBS39_P and ID1 reference profiles and sample fits](../analysis/release-comparison/extended-contexts/README.md). |

The figure pipeline starts from archived fits for the published attribution benchmarks and PCAWG decompositions.
Their linked runners regenerate those numerical inputs.
The driver and motif stages repeat their statistical analyses during the paper build.

## Recreate the compact inputs

The [preparation script](analysis/prepare_inputs.py) reads the original BGI MAF, hg19 genome, and 99 tumor spectra.
The [BGI example](../analysis/release-comparison/cohort-examples/README.md) documents the study archive and SBS96 profile command.
Extract `blca_bgi/data_clinical_sample.txt` from the same archive.
The [extended-context example](../analysis/release-comparison/extended-contexts/README.md) prepares both published reference matrices.
Pass those files to the script from the repository root.

```sh
python manuscript/analysis/prepare_inputs.py bgi.maf hg19.2bit bgi-counts.tsv data_clinical_sample.txt reference-directory
```

The script verifies the MAF checksum and checks tumor counts per allele against the archived ranking.
The paper reuses the archived exact matches to Cancer Hotspots v2.
The [original comparison](../analysis/release-comparison/pcawg-penta/README.md) documents transcript matching and the hotspot workbook.
The [bladder study](https://doi.org/10.1038/ng.2798) supplies the clinical records.

## Interpret the new analyses

- [Rank recovery](figures/rank-recovery.tsv) reports recall, precision, and precision relative to frequency at equal budgets.
  Equal scores share the available positions. Two thousand paired gene bootstrap draws estimate uncertainty conditional on the scores.
- [Allele scores](figures/allele-scores.tsv.gz) give recurrence, empirical mutability, reconstructed mutability, and exact hotspot membership.
  The fitted spectrum preserves the total substitution count. Both backgrounds use MutaGene's exome trinucleotide opportunities.
- [Gene recovery](figures/rank-recovery.tsv) uses distinct carrier counts as its baseline and the bundled driver gene list as its reference.
  The exploratory score multiplies the smallest signature-based B-score by the number of observed missense changes and caps it at one.
- [Motif groups](figures/motif-stage.tsv) use the cohort median fraction of eligible substitutions matching the APOBEC motif (0.389).
  The higher group includes ties. The two-sided Fisher test compares T2–T4 with Ta/T1 across groups.
  The analysis uses one specimen per clinical identifier and does not adjust for mutation burden or sequencing design.
- [Motif controls](figures/motif-controls.tsv) use 200 independent null comparisons and 200 planted comparisons.
  Each arm contains 30 tumors with 100 eligible substitutions. Doubling the motif weight changes its fraction from 1/8 to 2/9.
  This tests the single-motif comparison. The PCAWG screen adjusts all 7,715 tests together.

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
