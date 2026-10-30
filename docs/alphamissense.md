# Add AlphaMissense evidence

AlphaMissense annotation adds functional evidence without changing MutaGene's B-scores, q-values, or driver labels.
A pathogenicity score does not directly measure cancer-driver activity.

## Download and index a reference

Install the optional annotation dependency from this source checkout:

```sh
pip install '.[annotations]'
mutagene fetch alphamissense -g hg38
```

The command stores the canonical table under `~/.mutagene/alphamissense/hg38/`.
Use `-g hg19` for the other supported assembly.
Each compressed table requires approximately 600–650 MB.
The command resumes partial downloads, verifies size and MD5, and creates a tabix index.
It writes a provenance manifest with SHA256 hashes, assembly, source, rights, and verification time.
Repeated fetches verify the cached table before reuse.

The pinned source is [Zenodo record 10813168](https://zenodo.org/records/10813168), with CC-BY-4.0 terms.
The table accompanies [Cheng et al. (2023)](https://doi.org/10.1126/science.adg7492).
The canonical transcripts use GENCODE V27 for hg19 and V32 for hg38.
The cache does not contain model weights or perform new AlphaMissense inference.

## Annotate a ranking

The genome assembly must match the mutation coordinates and AlphaMissense table.
The MAF must include genomic alleles, transcript identifiers, and protein changes.
Run:

```sh
mutagene rank -i cohort.maf.gz -g hg38 \
  --alphamissense "$HOME/.mutagene/alphamissense/hg38/AlphaMissense_hg38.tsv.gz" \
  --alphamissense-genome hg38 -o ranked.tsv
```

The command adds `alphamissense_score`, `alphamissense_class`, and `alphamissense_status` columns.
Matching rows also include the source transcript and UniProt accession.
Versioned transcript identifiers require exact matches.
Unversioned matches carry `matched_unversioned_transcript` status.
Transcript or protein mismatches receive no score.
Ambiguous genomic mappings receive no score.
Missing coordinates receive `missing_coordinates` status.
Missing predictions remain missing; the implementation never converts them to zero.

The reader also accepts a separately prepared tabix-indexed table with the same columns.
Its records must identify the requested assembly.
Only the fetch workflow provides the pinned source checksum guarantee.
Preserve the source manifest when sharing a cache.

## Query without loading the full table

```python
from mutagene.io.alphamissense import AlphaMissenseLookup

with AlphaMissenseLookup("AlphaMissense_hg38.tsv.gz", "hg38") as lookup:
    result = lookup.lookup("1", 69094, "G", "T", "ENST00000335137.4", "V2L")
```

Genomic positions use one-based coordinates.
The example expects score `0.2937` and class `likely_benign`.
That example validates lookup behavior, not clinical interpretation.

## Evaluate a future combined model

A combined model requires cancer-specific labels and an independent evaluation plan.
Use these comparisons before changing any ranking:

1. Freeze MutaGene-only, AlphaMissense-only, and combined-model definitions.
2. Define driver and passenger labels independently of the candidate predictors.
3. Separate genes between training, calibration, and final evaluation.
4. Stratify by cancer type and oncogene or tumor-suppressor role.
5. Audit overlap with predictor training data and label-source circularity.
6. Compare precision-recall performance, recall at fixed false-discovery rates, and calibration.
7. Report annotation coverage and missingness for every comparison.
8. Assess temporal or independent-cohort transfer before prospective claims.

A calibrated supervised model can combine features while preserving separate evidence columns.
Do not multiply a p-value or B-score by a pathogenicity score.
The [2024 cancer assessment](https://doi.org/10.1016/j.jmb.2024.168644) explains the transfer problem.
[The 2026 MutaPheno study](https://doi.org/10.1073/pnas.2524289123) motivates independent evaluation with unseen proteins.
This implementation does not claim a measured driver-prediction gain.
