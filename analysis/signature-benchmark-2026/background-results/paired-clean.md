### 100 mutations

| Method | Total attribution error | Thresholded known error | Background mass | Active known | Cosine |
|---|---:|---:|---:|---:|---:|
| SigProfilerAssignment | 0.33306 | 0.32290 | 0.00000 | 3.06 | 0.87616 |
| no-background-BIC | 0.49488 | 0.46367 | 0.00000 | 3.55 | 0.88164 |
| block-BIC | 0.49488 | 0.46367 | 0.00000 | 3.55 | 0.88164 |
| joint-BIC | 0.53886 | 0.48497 | 0.03476 | 3.39 | 0.88343 |
| fixed-BIC | 0.55081 | 0.40169 | 0.27134 | 2.33 | 0.87847 |
| MLE | 0.67419 | 0.44202 | 0.05488 | 2.50 | 0.90830 |

### 2,000 mutations

| Method | Total attribution error | Thresholded known error | Background mass | Active known | Cosine |
|---|---:|---:|---:|---:|---:|
| joint-BIC | 0.04524 | 0.04512 | 0.00022 | 4.22 | 0.99084 |
| block-BIC | 0.04649 | 0.04649 | 0.00000 | 4.25 | 0.99084 |
| no-background-BIC | 0.04649 | 0.04649 | 0.00000 | 4.25 | 0.99084 |
| SigProfilerAssignment | 0.05121 | 0.05115 | 0.00000 | 3.70 | 0.98972 |
| fixed-BIC | 0.08339 | 0.06057 | 0.04565 | 3.94 | 0.99107 |
| MLE | 0.23824 | 0.21778 | 0.02083 | 14.98 | 0.99340 |

### 50,000 mutations

| Method | Total attribution error | Thresholded known error | Background mass | Active known | Cosine |
|---|---:|---:|---:|---:|---:|
| no-background-BIC | 0.00661 | 0.00661 | 0.00000 | 4.33 | 0.99964 |
| block-BIC | 0.00661 | 0.00661 | 0.00000 | 4.33 | 0.99964 |
| joint-BIC | 0.00706 | 0.00687 | 0.00039 | 4.33 | 0.99964 |
| SigProfilerAssignment | 0.00756 | 0.00756 | 0.00000 | 4.30 | 0.99963 |
| fixed-BIC | 0.01365 | 0.00931 | 0.00867 | 4.19 | 0.99965 |
| MLE | 0.05031 | 0.04788 | 0.00463 | 22.70 | 0.99974 |
