# PCAWG signatures and motifs at pentanucleotide resolution

The same 146 high-burden PCAWG tumors contribute 28,088,341 substitutions across 1,536 channels.
Collapsing the outer bases reproduces all 14,016 entries in the original SBS96 matrix.
The [exposure plot](pcawg-exposures.pdf) shows individual sparse BIC fits against 60 SignatureAnalyzer references.
The [exposures](exposures.tsv) identify each tumor and signature.
The sparse fits assign SBS3_P to 11 tumors. The largest exposure is 5.41% in Stomach-AdenoCA::SP105018.
Full-catalog maximum likelihood gives 5.32% in the same tumor.
Every fit passes the numerical optimality check for its method.
[Full-catalog exposures](mle-exposures.tsv.gz) · [All 292 fit checks](fit-checks.tsv)

## Compare a candidate with a supplied catalog

The [spectra](candidate-comparison.pdf) compare an extra component with its nearest supplied reference and the best nonnegative mixture.
The full-catalog candidate most closely resembles SBS17b_P (cosine 0.953).
The fitted mixture matches at 0.962 and assigns 70.5% to SBS17b_P and 17.5% to SBS28_P.
The mean reduction in held-out error is 0.000035. Its interval includes zero.

Removing SBS7a_S tests recovery of a known omission.
The resulting component matches SBS7a_S at 0.999. Its nearest remaining reference matches at 0.847.
Mean held-out error falls by 0.0630 (95% interval 0.0421 to 0.0838).
This control demonstrates recovery of a missing signature. The full-catalog experiment provides little evidence for an additional biological process.

[Candidate spectra](candidate-profiles.tsv) · [Every reference match and mixture weight](candidate-matches.tsv) · [Held-out tumors](held-out.tsv) · [Summary](candidate-summary.tsv)

We fix the supplied profiles and learn one extra nonnegative component by alternating least squares.
Each training tumor contributes its channel proportions.
Positive residuals from the reference-only fit initialize the component.
The split reserves one-third of each histology, rounded down, for testing (104 training tumors, 42 test tumors).
A binomial split assigns test mutations to exposure fitting or evaluation with equal probability.
Both models use nonnegative least squares for exposures and normalize predictions to unit mass.
Error equals half the sum of absolute differences between predicted and observed channel proportions.
The improved-tumor count excludes error reductions below numerical precision (at most 1e-12).
The interval uses 5,000 paired tumor bootstrap draws and conditions on the learned component.
The supplied reference catalog also uses PCAWG data. Independent biological validation requires a separate cohort.

For comparison, unrestricted rank-five NMF learns five profiles without the catalog.
Ten random starts compete on squared reconstruction error.
The nearest matches range from 0.867 to 0.992. [All components](extracted-signatures.tsv) · [Five best matches](signature-matches.tsv)

## Identify enriched motifs

The screen tests seven curated motifs on both strands and all 1,536 exact substitution channels.
Each of five histologies has at least ten tumors. Every comparison uses all remaining tumors as its reference.
Each tumor contributes one motif fraction within compatible substitution classes.
A one-sided Mann–Whitney test compares tumor fractions. Benjamini–Hochberg adjusts all 7,715 tests together.
Fold enrichment compares the two medians. Genomic opportunity can also affect motif frequencies.
The site-based MutaGene motif command uses local genomic opportunities and answers a different question.

Melanoma shows 8.00-fold enrichment for CT[T>A]TC (65 versus 81 tumors).
Colorectal tumors show 4.98-fold enrichment for [C>T]G (16 versus 130 tumors).
The [full screen](motifs.tsv) contains motif identities, counts of tumors, group medians, fold changes, and adjusted p-values.

## Rank drivers and compare known hotspots

The BGI bladder MAF contains 99 tumors. Ranking uses each observed protein change and its background mutability.
The cohort supplies its own background and recurrence.
Tumor B10 is the first sample by lexical order. Its ranking uses the 240-tumor Bladder_Urothelial_Carcinoma background.
An SBS count matrix lacks the protein identities that ranking requires.

| Cohort comparison | Driver group | Other scored group | Fold | One-sided Fisher p |
|---|---:|---:|---:|---:|
| Cancer Hotspots v2 exact missense changes | 57/2,451 | 30/4,188 | 3.25 | 4.87 × 10⁻⁸ |
| Bundled common driver genes | 19/2,074 | 15/2,673 | 1.63 | 0.103 |

Each denominator counts distinct missense changes or distinct genes.
A gene enters the Driver group if any scored missense change meets the Driver threshold.
Exact matches require the gene, reference residue, position, and alternate residue.
Reference transcript numbering affects matches.
B10 has 45 Driver and ten other scored missense changes. Neither group contains an exact hotspot match.
Both methods use recurrence. Study overlap with the hotspot catalog and reference cohort remains unknown.
The comparison measures catalog agreement. Unlisted changes have unknown driver status.

[Cohort ranking](bgi-ranked.tsv.gz) · [B10 ranking](B10-ranked.tsv) · [Cohort enrichment](cohort-enrichment.tsv) · [B10 enrichment](sample-enrichment.tsv) · [Figure](motifs-and-drivers.pdf)

## Reproduce the signature and motif analyses

Use the benchmark Python environment, MutaGene, and R from the repository root.
The [prepared count matrix](pcawg-counts.tsv.gz) contains integer SBS1536 counts.
The [upstream PCAWG matrix](https://github.com/getzlab/SignatureAnalyzer/blob/f205be05f6e180bd04c3beaffc354b73ae3ccf9f/PCAWG/INPUT_SignatureAnalyzer/lego1536.PAN.SNV.091217.RData) has SHA-256 `f473f3787b07d1d4da9decc8d353d9a56d94f0752b056103b6e3f6c0ffb8ea5c`.
To regenerate the count matrix, prepare the [original SBS96 counts](../cohort-examples/README.md) first.
Run `Rscript prepare.R source.RData original-SBS96.tsv counts.tsv` from this directory.
The script matches specimen identifiers and checks every collapsed count.

The [extended-context example](../extended-contexts/README.md) prepares the 60-signature reference.
These commands expect that prepared reference in the working directory.

```sh
d=analysis/release-comparison/pcawg-penta
ref=SBS1536-reference.tsv
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
gzip -dc "$d/pcawg-counts.tsv.gz" > counts.tsv
mutagene signature --context SBS1536 -f MATRIX -i counts.tsv --reference-matrix "$ref" \
    -o "$d/exposures.tsv" --diagnostics-out fit-diagnostics.json
mutagene discover -i counts.tsv --reference-matrix "$ref" -o "$d" --withhold SBS7a_S
python "$d/../cohort-examples/extract.py" counts.tsv --reference "$ref" --output "$d"
python "$d/motifs.py" counts.tsv "$d/motifs.tsv"
python "$d/figure.py" "$ref"
```

To reproduce the full-catalog exposures, repeat the signature command with `-m MLE` and separate output paths.
For a custom catalog, replace the reference path and omit the optional control.
The candidate comparison accepts labeled nonnegative reference profiles and integer tumor counts with matching channels.
The command creates the output directory and writes the learned profile and every reference match.

## Reproduce the driver comparison

Prepare the [BGI MAF](../cohort-examples/README.md) and an hg19 genome.
The single-tumor example also requires the MutaGene cohort archive.
Download the [Cancer Hotspots v2 workbook](https://www.cancerhotspots.org/files/hotspots_v2.xls).
Its SHA-256 is `7856cb650b7f2d4f82448c2e08c275eb5f828adcfcdd13f33cc6f749288014d9`.
The converter requires pandas and xlrd. It retains 2,464 distinct missense changes on the SNV-hotspots sheet.

```sh
d=analysis/release-comparison/pcawg-penta
python "$d/hotspots.py" hotspots_v2.xls hotspots.tsv
mutagene rank -i bgi.maf -g hg19 -o bgi-ranked.tsv --known-variants hotspots.tsv \
    --enrichment-out "$d/cohort-enrichment.tsv"
python -c 'import pandas as pd; t=pd.read_csv("bgi.maf", sep="\t", comment="#", dtype=str, keep_default_na=False); t[t.Tumor_Sample_Barcode=="B10"].to_csv("B10.maf",sep="\t",index=False)'
mutagene rank -i B10.maf -g hg19 -c Bladder_Urothelial_Carcinoma --cohorts-file cohorts.tar.gz \
    -o B10-ranked.tsv --known-variants hotspots.tsv --enrichment-out "$d/sample-enrichment.tsv"
```
