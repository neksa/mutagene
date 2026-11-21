# Cohort decomposition and treatment outcomes

The [figure](cohorts.pdf) shows MutaGene exposures for PCAWG tumors and published SigMA groups for two independent treatment studies.
Each panel uses a different cohort.

## PCAWG exposures

Medo et al. provide SBS96 counts for 146 PCAWG whole-genome catalogs with at least 50,000 substitutions.
The counts and histology labels come from the [pinned SigFitTest input](https://github.com/8medom/SigFitTest/blob/e94bac6e9b04fb13500d1904b5e0e8d58b00d8ad/real%20mutational%20catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat).
MutaGene fits each tumor with joint sparse BIC, the [67-signature catalog](../performance/reference.tsv), and six substitution-specific background components.
The [exposure table](pcawg-weights.tsv) reports each reference weight and total fitted background mass.
The stacked bars show SBS3 separately and highlight five other signatures by mean exposure.
Each tumor contributes equally to that ranking.
MutaGene assigns nonzero SBS3 weights to 19 of 146 tumors.
Stomach-AdenoCA::SP105018 has the largest fitted weight (17.6%).
A zero weight means that the sparse fit assigns zero exposure.
The PCAWG analysis uses continuous exposures. The treatment studies use their own SigMA classifications.
Within each group, the figure sorts tumors by decreasing exposure to the five leading signatures, in rank order.
Lung combines adenocarcinoma and squamous cell carcinoma. Other includes the remaining histologies.
The burden threshold limits generalization to tumors with fewer substitutions.

We also check survival coverage in the cBioPortal PCAWG clinical tables.
Joining the 146 specimen identifiers to SAMPLE_ID identifies 131 distinct patients.
Only six have both OS_MONTHS and OS_STATUS.
They span three histologies. This subset offers little information for survival comparisons by signature exposure.
[Sample records](https://github.com/cBioPortal/datahub/blob/2a86a1f9289523cb5bdcceec5a2c61cf6a93ca26/public/pancan_pcawg_2020/data_clinical_sample.txt) · [Patient records](https://github.com/cBioPortal/datahub/blob/e7f71c358b4f7ecf9fc04f29823c2a487b4d32e1/public/pancan_pcawg_2020/data_clinical_patient.txt)

## Extract candidate SBS96 signatures

The [extraction script](extract.py) uses scikit-learn 1.9.1 and the same 146 PCAWG profiles as panel A.
Each tumor contributes its SBS96 proportions. The reference catalog enters only after extraction.
We fix the rank at five for this example. The script uses standard NMF with squared Frobenius loss.
It selects the fit with the smallest reconstruction error across ten random starts (seeds 0 through 9).
It scales each [component profile](extracted-signatures.tsv) to unit mass and orders components by mean activity.
The [matches](signature-matches.tsv) give the largest cosine similarity across the 67 reference signatures for each component.
A weak match alone cannot establish a new mutational process.
The original study uses nonsmooth NMF to derive MutaGene-5 and MutaGene-10. [Original study](https://doi.org/10.1093/nar/gkx367)

After preparing the PCAWG counts below, run the following command from the repository root.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python analysis/release-comparison/cohort-examples/extract.py /tmp/pcawg-counts.tsv
```

## Test a motif in the example tumor

The hg19 assembly must be available locally.
Run this command from the repository root to reproduce the curated motif screen.

```sh
mutagene motif -i sample1.maf -g hg19 --strand A -o motifs.tsv
```

The CpG motif `[C>T]G` has enrichment 2.11727 and adjusted q-value 0.0118927 in sample TCGA-50-6593-01A-11D-1753-08.
The motif window defaults to 50 bases on each side.

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
These percentages correspond to 14 and four responders, respectively. The plot reads the [group counts](petremac.tsv).
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

## Compare the five components in independent cohorts

We repeat the same rank-five extraction on 99 bladder tumors from the BGI study.
The source MAF contains 10,655 single-base substitutions.
The median tumor contributes 79 substitutions (range 4 to 464).
This single-organ collection differs from the high-burden PCAWG discovery panel.
The coding calls and sequencing designs also affect comparisons. [BGI study](https://doi.org/10.1038/ng.2798)

We compare each PCAWG component with all five [BGI components](bgi-signatures.tsv).
We also compare each component with 315 published Genomics England (GEL) profiles and 249 published Hartwig profiles.
These organ-specific catalogs include profiles from independent extractions. [GEL and Hartwig study](https://doi.org/10.1126/science.abl9283)
The [results](replication.tsv) give the highest cosine similarity within each collection.
Larger collections provide more chances for close matches. These maxima describe spectral agreement.
The [comparison figure](replication.pdf) shows the four spectra for component 3 with identical axes.

Component 3 matches a BGI component at 0.797, GEL SBS2+13 at 0.940, and Hartwig SBS13 at 0.929.
A two-signature NNLS fit uses the unit-mass COSMIC v3 SBS2 and SBS13 profiles.
The normalized weights are 44.4% SBS2 and 55.6% SBS13. The mixture gives a cosine similarity of 0.908.
These results support a combination of known patterns. They provide no independent evidence for a new biological process.

To repeat the BGI extraction, use the benchmark environment, hg19, and the [study archive](https://datahub.assets.cbioportal.org/blca_bgi.tar.gz).
The source MAF has SHA-256 `c09b911056ba5d5f750df589a856706a1eb782290b78bdf32b68521fcfffe7d7`.
The comparison script retrieves version 2.03 of both external catalogs from a pinned commit.

```sh
examples=analysis/release-comparison/cohort-examples
mkdir -p /tmp/bgi-extraction
tar -xzOf blca_bgi.tar.gz blca_bgi/data_mutations.txt > /tmp/bgi.maf
mutagene profile --per-sample -f MAF -g hg19 -i /tmp/bgi.maf -o /tmp/bgi-counts.tsv
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$examples/extract.py" /tmp/bgi-counts.tsv --output /tmp/bgi-extraction
python "$examples/replicate.py" /tmp/bgi-extraction/extracted-signatures.tsv
```
