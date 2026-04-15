# Endometriosis Gene Signature Discovery

## Overview

This project implements a complete and reproducible pipeline for identifying a compact gene expression signature associated with endometriosis using publicly available transcriptomic data.

The workflow combines data preprocessing, statistical filtering, machine learning, feature stability analysis, and biological annotation to extract a minimal set of genes that collectively capture disease-related signal.

The goal is not only predictive performance, but also interpretability and structural understanding of the underlying transcriptional patterns.

---

## Dataset

- Source: Gene Expression Omnibus (GEO)
- Study: GSE51981
- Platform: Affymetrix GPL570 microarray

The dataset consists of endometrial tissue samples from patients with and without endometriosis, including metadata describing tissue type and disease severity.

Files needed:

GPL570-tbl-1.txt

GSE51981_family.xml

GSE51981_series_matrix.txt

---

## Methodology

The pipeline is structured in sequential stages:

### 1. Data Preparation

- Gene expression matrix is loaded and transposed (samples × genes)
- Metadata is extracted from GEO XML (MINiML format)
- Samples are labeled:
  - `1` → Endometriosis
  - `0` → Non-endometriosis
- Expression data and metadata are merged into a unified dataset

### 2. Cohort Filtering

To reduce confounding, a clean subset is constructed:

- Included:
  - Endometriosis (all severities)
  - Healthy controls without uterine pelvic pathology
- Excluded:
  - Controls with other pathologies

This step significantly improves signal clarity.

---

### 3. Baseline Model

A logistic regression model is used with:

- Variance threshold filtering
- Standardization
- 5-fold cross-validation

This provides a robust baseline for predictive performance.

---

### 4. Feature Stability Analysis

To identify meaningful predictors:

- Model is trained across multiple CV folds
- Top features (by coefficient magnitude) are extracted per fold
- Feature occurrence is counted across folds

This produces a **stability score** for each probe, identifying consistently selected features.

---

### 5. Probe Annotation

Affymetrix probe IDs are mapped to gene symbols using the GPL570 annotation table.

This step bridges machine learning outputs to biological interpretation.

---

### 6. Feature Ranking

Features are ranked based on:

- Stability (frequency across folds)
- Model importance (coefficient magnitude)

This produces a prioritized list of candidate genes.

---

### 7. Signature Extraction

A gene signature is defined as the top-ranked subset of probes.

The signature is evaluated independently using cross-validation to assess whether a reduced feature set retains predictive power.

---

### 8. Signature Size Optimization

Performance is evaluated as a function of signature size:

- Signature size is varied (2–14 genes)
- ROC-AUC is computed for each size
- The minimal size achieving maximal performance is selected

---

## Key Results

- Full model (high-dimensional):
  - ROC-AUC ≈ 0.93

- Compact signature (~7 genes):
  - ROC-AUC ≈ 0.95

- Performance plateaus beyond ~7 genes

### Interpretation

- The disease signal is **low-dimensional**
- A small subset of genes captures most of the predictive information
- Additional features are largely redundant

---

## Identified Signature (Dataset-Specific)

The final signature consists of a small set of probes mapping to genes including:

- CTU2
- ZNF24
- NT5DC3
- LINC00925
- HMGN3-AS1
- ZNF568
- C11orf54

These genes form a **distributed transcriptional signature**, rather than a single dominant biomarker.

---

## Project Structure

```
.
├── data/
│   ├── raw/                # GEO raw data and platform annotation
│   └── processed/          # cleaned datasets and results
├── scripts/
│   ├── build_dataset.py    # data preprocessing pipeline
│   └── run_signature.py    # signature discovery pipeline
├── src/
│   └── endometriosis_signature/
│       ├── dataset.py      # data loading and merging
│       ├── modeling.py     # ML pipeline and evaluation
│       ├── annotations.py  # probe-to-gene mapping
│       └── constants.py    # configuration values
└── README.md
```

---

## Usage

Install the project:

```
pip install -e .
```

Run dataset construction:

```
python scripts/build_dataset.py
```

Run signature discovery:

```
python scripts/run_signature.py
```

---

## Outputs

Generated files include:

- `endometriosis_clean_filtered.parquet`
- `stable_probe_counts.csv`
- `stable_mapped_probes.csv`
- `ranked_signature_probes.csv`
- `signature_size_sweep.csv`

These files capture intermediate and final results of the pipeline.

---

## Limitations

- Results are dataset-specific and not externally validated
- Microarray probe annotation is incomplete for some probes
- Cohort size is limited, which may affect generalizability

---

## Future Work

- External validation on independent datasets
- Pathway enrichment analysis
- Sparse models (e.g., L1 regularization)
- Integration with RNA-seq data
- Biological interpretation of identified genes

---

## Conclusion

This project demonstrates that endometriosis-associated transcriptional patterns can be captured by a compact and stable gene signature. The approach highlights the importance of combining machine learning with stability analysis and biological annotation to extract meaningful structure from high-dimensional data.
