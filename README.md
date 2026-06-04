# Endometriosis Gene Signature Discovery

This repository contains a reproducible transcriptomics workflow for identifying a compact gene-expression signature associated with endometriosis using public GEO datasets.

The project is designed as a complete bioinformatics pipeline: it starts from raw expression data and sample metadata, builds a cleaned case-control dataset, trains interpretable logistic-regression models, ranks stable probes, evaluates signature size, and performs external validation across independent datasets.

## Project aim

The main goal is not only to classify endometriosis and control samples, but to understand how much of the disease-associated transcriptomic signal can be captured by a small and interpretable set of probes.

The project focuses on three questions:

1. Can endometriosis and control samples be separated using public transcriptomic data?
2. Can this signal be reduced to a compact and stable probe signature?
3. How well does the signature transfer to external datasets generated on different platforms?

## Repository structure

```text
endometriosis-signature/
├── data/
│   ├── raw/
│   │   └── raw GEO files supplied by the user
│   └── processed/
│       ├── cleaned datasets
│       ├── ranked probe tables
│       ├── cross-validation metrics
│       └── signature-size sweeps
├── results/
│   └── figures/
│       ├── ROC curves
│       ├── signature-size curves
│       └── exploratory stable-probe sweeps
├── scripts/
│   ├── build_dataset.py
│   ├── run_signature.py
│   ├── run_signature_l1.py
│   ├── plot_signature_size.py
│   ├── plot_stable_probe_sweep.py
│   ├── validate_gse135485.py
│   └── validate_gse25628.py
├── src/
│   └── endometriosis_signature/
│       ├── annotation.py
│       ├── constants.py
│       ├── dataset.py
│       └── modeling.py
├── pyproject.toml
└── README.md
```

## Data

The main training dataset is based on GSE51981, generated on the GPL570 Affymetrix Human Genome U133 Plus 2.0 Array platform.

External validation is performed on independent datasets, including GSE135485 and GSE25628.

Raw GEO files are expected in:

```text
data/raw/
```

They are not included in the repository.

## Dataset construction

The dataset-building step is handled by:

```bash
python scripts/build_dataset.py
```

This script performs the following steps:

1. Loads the GSE51981 expression matrix.
2. Transposes the matrix so that samples are rows and probes are columns.
3. Parses sample-level metadata.
4. Builds a binary target variable for endometriosis status.
5. Merges expression values and metadata.
6. Filters the cohort to a cleaner case-control comparison.
7. Saves the processed dataset as a Parquet file.

In the current run, the script produced:

```text
Full dataset: 148 samples, 54679 columns
Clean dataset: 109 samples, 54679 columns
Clean labels: 75 endometriosis, 34 controls
```

## Modeling strategy

The main model is an interpretable logistic-regression pipeline implemented with scikit-learn.

The pipeline contains three steps:

1. Variance filtering
2. Standardization
3. Logistic regression

Low-variance probes are removed because probes that are nearly constant across all samples cannot contribute meaningfully to case-control separation. This step is unsupervised and does not use the class labels.

Standardization is then applied so that all probes are placed on a comparable scale before fitting the model.

Logistic regression is used because it provides a transparent linear classifier whose coefficients can be inspected and ranked. This is useful for transcriptomic signature discovery, where interpretability is as important as classification performance.

## Cross-validation

Model performance is evaluated using 5-fold cross-validation with ROC-AUC as the main metric.

In each split:

1. The model is trained on 80 percent of the samples.
2. It is evaluated on the remaining 20 percent.
3. ROC-AUC is computed on the held-out fold.

The baseline L2 logistic-regression model reached:

```text
Baseline AUC: 0.933 ± 0.058
```

## Stable probe ranking

Feature stability is estimated across cross-validation folds.

For each fold:

1. The model is trained on the training split.
2. Logistic-regression coefficients are extracted.
3. Probes with the strongest coefficients are recorded.
4. Probe recurrence across folds is counted.

Probes that repeatedly appear among the strongest coefficients are considered more stable candidates.

The stable probes are mapped to gene symbols using GPL570 annotation and ranked by:

1. Cross-validation recurrence count
2. Absolute logistic-regression coefficient magnitude

The final mapped ranked signature contains 7 probes:

```text
CTU2
ZNF24
LINC00925
NT5DC3
HMGN3-AS1
C11orf54
ZNF568
```

The corresponding probe-level table is saved as:

```text
data/processed/ranked_signature_probes.csv
```

## Signature-size analysis

The script:

```bash
python scripts/run_signature.py
```

also evaluates how performance changes as ranked probes are progressively added to the model.

The corrected signature-size sweep uses only real available ranked probes. Since the final mapped ranked signature contains 7 probes, the final curve tests signature sizes from 2 to 7 probes.

Current results:

```text
2 probes -> AUC: 0.899 ± 0.109
3 probes -> AUC: 0.899 ± 0.109
4 probes -> AUC: 0.899 ± 0.109
5 probes -> AUC: 0.920 ± 0.079
6 probes -> AUC: 0.933 ± 0.065
7 probes -> AUC: 0.948 ± 0.047
```

The best observed performance is obtained with the full 7-probe mapped signature.

Outputs:

```text
data/processed/signature_size_sweep.csv
results/figures/signature_size_curve.png
```

## Exploratory stable-probe sweep

A second exploratory analysis is included:

```bash
python scripts/plot_stable_probe_sweep.py
```

This analysis uses the broader pool of 64 stable probes from:

```text
data/processed/stable_probe_counts.csv
```

The purpose is to inspect whether the broader stable-probe pool contains additional predictive signal and whether performance appears to saturate as more probes are added.

This analysis should be interpreted carefully. The stable-probe pool was already selected using the dataset, so the resulting curve is optimistic and should not be treated as an unbiased validation estimate. It is included as a diagnostic ranking analysis, not as the main performance claim.

Outputs:

```text
data/processed/stable_probe_size_sweep.csv
results/figures/stable_probe_size_curve.png
```

## L1 logistic-regression analysis

The repository also includes an L1-regularized logistic-regression analysis:

```bash
python scripts/run_signature_l1.py
```

The L1 model acts as a sparse embedded feature selector. While the L2 model distributes weight across correlated transcriptomic features, the L1 model pushes many coefficients exactly to zero and identifies a smaller core signature.

This analysis is used to ask a complementary question:

```text
What is the smallest set of probes that still captures a large fraction of the signal?
```

The L1 sweep evaluates different regularization strengths and saves the selected probes and genes for each value of C.

Output:

```text
data/processed/signature_l1_sweep.csv
```

## External validation

External validation scripts test whether the discovered signature transfers to independent datasets.

```bash
python scripts/validate_gse135485.py
python scripts/validate_gse25628.py
```

These validations are important because transcriptomic signatures can perform well within one dataset but lose performance across studies, technologies, or preprocessing pipelines.

The external validation results are therefore interpreted as a test of transferability rather than as a direct replacement for the internal cross-validation results.

## How to run the project

Create and activate an environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .
```

Place the required GEO files in:

```text
data/raw/
```

Then run:

```bash
python scripts/build_dataset.py
python scripts/run_signature.py
python scripts/plot_signature_size.py
python scripts/run_signature_l1.py
python scripts/plot_stable_probe_sweep.py
python scripts/validate_gse135485.py
python scripts/validate_gse25628.py
```

## Main outputs

```text
data/processed/endometriosis_clean_filtered.parquet
data/processed/stable_probe_counts.csv
data/processed/stable_mapped_probes.csv
data/processed/ranked_signature_probes.csv
data/processed/signature_size_sweep.csv
data/processed/signature_l1_sweep.csv
data/processed/stable_probe_size_sweep.csv
results/figures/signature_size_curve.png
results/figures/stable_probe_size_curve.png
```

## Current interpretation

The cleaned GSE51981 cohort contains a strong transcriptomic signal separating endometriosis and control samples.

A compact 7-probe mapped signature achieved the best internal 5-fold cross-validation performance among the tested mapped signature sizes.

The broader exploratory stable-probe sweep suggests that the larger stable-probe pool contains additional disease-separation signal, but that analysis is optimistic because the candidate pool was defined before the sweep. The main reported signature-size result is therefore the corrected 7-probe mapped-signature curve.

## Limitations

This project is intended as a reproducible exploratory transcriptomics workflow, not as a clinical diagnostic model.

Important limitations include:

- Small sample size after cohort filtering
- Class imbalance between endometriosis and controls
- Potential batch and platform effects
- Probe-level rather than fully gene-collapsed modeling
- Optimism in analyses where feature selection is performed before cross-validation
- Need for larger independent validation cohorts

A stronger future version would implement nested cross-validation, where probe ranking and feature selection are repeated inside each training fold before evaluation on the held-out fold.

## Technical summary

This project demonstrates:

- GEO transcriptomic data processing
- Microarray metadata parsing
- Cohort cleaning
- Probe annotation
- Scikit-learn modeling pipelines
- Variance filtering
- Standardization
- L2 and L1 logistic regression
- Cross-validation with ROC-AUC
- Coefficient-based feature ranking
- Stability-based probe selection
- Signature-size evaluation
- External dataset validation
- Reproducible Python project organization
