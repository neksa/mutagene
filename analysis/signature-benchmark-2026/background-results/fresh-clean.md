### 100 mutations

| Method | Total attribution error | Thresholded known error | Background mass | Active known | Cosine |
|---|---:|---:|---:|---:|---:|
| SigProfilerAssignment | 0.28875 | 0.28414 | 0.00000 | 2.84 | 0.86604 |
| block-BIC | 0.46893 | 0.43835 | 0.00000 | 3.05 | 0.86788 |
| no-background-BIC | 0.46893 | 0.43835 | 0.00000 | 3.05 | 0.86788 |
| joint-BIC | 0.50179 | 0.45346 | 0.03829 | 3.16 | 0.87217 |
| fixed-BIC | 0.52158 | 0.36873 | 0.28504 | 2.12 | 0.86638 |
| MLE | 0.64675 | 0.42486 | 0.06651 | 2.45 | 0.89806 |

### 2,000 mutations

| Method | Total attribution error | Thresholded known error | Background mass | Active known | Cosine |
|---|---:|---:|---:|---:|---:|
| no-background-BIC | 0.05537 | 0.05537 | 0.00000 | 4.50 | 0.99034 |
| block-BIC | 0.05537 | 0.05537 | 0.00000 | 4.50 | 0.99034 |
| SigProfilerAssignment | 0.05628 | 0.05627 | 0.00000 | 3.75 | 0.98873 |
| joint-BIC | 0.05823 | 0.05739 | 0.00168 | 4.50 | 0.99035 |
| fixed-BIC | 0.09446 | 0.06721 | 0.05444 | 4.05 | 0.99046 |
| MLE | 0.25692 | 0.23568 | 0.02111 | 15.22 | 0.99327 |

### 50,000 mutations

| Method | Total attribution error | Thresholded known error | Background mass | Active known | Cosine |
|---|---:|---:|---:|---:|---:|
| joint-BIC | 0.00591 | 0.00591 | 0.00000 | 4.34 | 0.99960 |
| no-background-BIC | 0.00604 | 0.00604 | 0.00000 | 4.39 | 0.99960 |
| block-BIC | 0.00604 | 0.00604 | 0.00000 | 4.39 | 0.99960 |
| SigProfilerAssignment | 0.00747 | 0.00747 | 0.00000 | 4.36 | 0.99959 |
| fixed-BIC | 0.01378 | 0.00939 | 0.00878 | 4.30 | 0.99961 |
| MLE | 0.05577 | 0.05304 | 0.00521 | 23.55 | 0.99972 |
