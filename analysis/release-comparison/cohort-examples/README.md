# Cohort decomposition and treatment outcomes

The [figure](cohorts.pdf) shows MutaGene exposures for PCAWG tumors and published SigMA groups for two independent treatment studies.
Each panel uses a different cohort.

## PCAWG exposures

Medo et al. provide SBS96 counts for 146 PCAWG whole-genome catalogs with at least 50,000 substitutions.
The counts and histology labels come from the [pinned SigFitTest input](https://github.com/8medom/SigFitTest/blob/e94bac6e9b04fb13500d1904b5e0e8d58b00d8ad/real%20mutational%20catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat).
MutaGene fits each tumor with joint sparse BIC, the [67-signature catalog](../performance/reference.tsv), and six substitution-specific background components.
The [exposure table](pcawg-weights.tsv) reports each reference weight and total fitted background mass.
Each tumor contributes equally to the ranking of signatures by mean exposure.
Within each group, the figure sorts tumors by decreasing exposure to the five leading signatures, in rank order.
Lung combines adenocarcinoma and squamous cell carcinoma. Other includes the remaining histologies.
The burden threshold limits generalization to tumors with fewer substitutions.

We also check survival coverage in the cBioPortal PCAWG clinical tables.
Joining the 146 specimen identifiers to SAMPLE_ID identifies 131 distinct patients.
Only six have both OS_MONTHS and OS_STATUS.
They span three histologies. This subset offers little information for survival comparisons by signature exposure.
[Sample records](https://github.com/cBioPortal/datahub/blob/2a86a1f9289523cb5bdcceec5a2c61cf6a93ca26/public/pancan_pcawg_2020/data_clinical_sample.txt) · [Patient records](https://github.com/cBioPortal/datahub/blob/e7f71c358b4f7ecf9fc04f29823c2a487b4d32e1/public/pancan_pcawg_2020/data_clinical_patient.txt)

## TOPACIO progression-free survival

The TOPACIO source workbook records 62 patients with recurrent ovarian cancer.
Its Figure 1 sheet includes SigMA status for 39 patients receiving niraparib plus pembrolizumab.
The [event counts](topacio.tsv) group observations by status, follow-up day, and censoring indicator.
A censor value of one denotes censoring. Zero denotes a progression-free survival event.
The Kaplan–Meier estimator includes all observations with known status. Crosses mark censoring.
The time axis uses 365.25/12 days per month.
The risk table counts patients under observation immediately before each indicated time.
R survival gives median progression-free survival of 152 days for positive patients and 66 days for negative patients.
The paper reports 5.0 and 2.2 months, respectively. [Study](https://doi.org/10.1038/s41467-020-15315-8) · [Source workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-15315-8/MediaObjects/41467_2020_15315_MOESM5_ESM.xlsx)

## PETREMAC objective response

Batalini et al. report objective response to neoadjuvant olaparib in 32 patients with primary triple-negative breast cancer.
Figure 4D reports 74% response among 19 panel-Sig3-positive patients and 31% among 13 panel-Sig3-negative patients.
These percentages correspond to 14 and four responders, respectively.
The figure shows 95% Wilson binomial intervals calculated from those counts.
Both treatment studies use SigMA to determine Sig3 status.
Neither study evaluates MutaGene predictions. [Study, cohort 2](https://doi.org/10.1158/1078-0432.CCR-22-0749)

## Reproduce the figure

Use the [benchmark environment](../performance/README.md) with openpyxl for the source workbook.
Run the commands from the repository root.
Download the linked PCAWG input and TOPACIO workbook to the indicated temporary paths.
The preparation script verifies both source checksums.

```sh
examples=analysis/release-comparison/cohort-examples
performance=analysis/release-comparison/performance
python "$examples/prepare.py" pcawg /tmp/pcawg.tsv /tmp/pcawg-counts.tsv
PYTHONPATH="$PWD" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$performance/fit.py" \
  MutaGene /tmp/pcawg-counts.tsv "$performance/reference.tsv" "$examples/pcawg-weights.tsv"
python "$examples/prepare.py" topacio /tmp/topacio.xlsx "$examples/topacio.tsv"
python "$examples/figure.py"
```

Use the last command alone to render the figure from the accompanying tables.
