# Scientific review of the application note

Claude reviewed revision `972665c1` against the manuscript, analysis code, and archived results.
agy reviewed the revised text, generated tables, and benchmark methods.

| Reviewer | Finding | Revision |
|---|---|---|
| Claude | The pooled attribution error conceals differences between background conditions. | The abstract qualifies the comparison. Page 8 gives both methods' errors under all four conditions. |
| agy | Methods account for unassigned signature weight differently. | Page 2 explains which methods receive zero background mass and how the benchmark scores deconstructSigs. |
| Claude | Gene bootstrap intervals omit variation between possible tie orders. | Page 7 reports exact recurrence tie ranges separately from uncertainty in tie-averaged precision. |
| agy | Two scoring models share the B-score name. | Page 8 identifies the empirical B-score as the source of the 57 hotspot matches. |
| Claude | High cosine similarity alone provides limited evidence about novelty. | Page 6 reports held-out gain and its conditioning assumptions for both discovery experiments. |
| agy | The external comparison mixes cohort extractions and reference catalogs. | Page 9 identifies the BGI extraction and the published GEL and Hartwig catalogs separately. |
| Claude | Clinical group sizes and the source of response counts need clarification. | Page 3 identifies the TOPACIO subset and the counts implied by PETREMAC percentages. |
| Claude | The headline benchmark lacks a pipeline dependency. | DVC now generates the pooled errors and condition table. It also tracks the PETREMAC counts. |

The revised driver table uses two sources of uncertainty.
A hypergeometric distribution describes the possible hotspot counts when recurrence ties receive random order within this cohort.
The paired gene bootstrap estimates precision differences across resampled genes and averages over ties within each sample.
These calculations answer different questions.
