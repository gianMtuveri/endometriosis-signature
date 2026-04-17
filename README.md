# Endometriosis Gene Signature Discovery

## Overview

This project builds a reproducible pipeline to identify a compact gene expression signature associated with endometriosis using public transcriptomic datasets.

The goal is not only to obtain good predictive performance, but to understand how much of the disease signal can be captured by a small and interpretable set of genes, and how well that signal transfers across independent datasets and experimental platforms.

---

## Datasets

### Training dataset
- GEO: GSE51981
- Platform: Affymetrix GPL570 microarray
- Task: endometriosis vs healthy controls

### External validation datasets
- GSE135485
  - Platform: RNA-seq
  - Purpose: cross-platform validation
- GSE25628
  - Platform: GPL571 microarray
  - Purpose: cross-study validation

---

## Pipeline

### 1. Data preparation
The training dataset is built by loading the GEO expression matrix, parsing metadata from XML, and merging sample-level annotations with gene expression values.

### 2. Cohort cleaning
The original study contains heterogeneous controls. To reduce confounding, the analysis keeps:
- patients with endometriosis
- healthy controls without uterine pelvic pathology

This produces a cleaner and more interpretable classification problem.

### 3. Baseline modeling
The main baseline uses:
- variance filtering
- standardization
- logistic regression
- 5-fold cross-validation

### 4. Feature stability analysis
Instead of trusting a single model fit, features are ranked according to how consistently they appear across cross-validation folds.

This identifies a stable subset of predictive probes rather than one-off correlations.

### 5. Annotation and ranking
Stable probes are mapped to gene symbols using GPL570 annotation, then ranked using:
- cross-validation stability
- coefficient magnitude

### 6. Signature extraction
A compact signature is derived from the top-ranked probes.

### 7. Signature size analysis
The number of genes is varied and model performance is tracked to determine where the signal saturates.

### 8. Sparse signature analysis
L1-regularized logistic regression is used to test whether the signature can be compressed further into an even smaller subset of genes.

### 9. External validation
The resulting signatures are tested on:
- an RNA-seq dataset
- an independent microarray dataset

This evaluates robustness across both studies and platforms.

---

## Main Results

### Summary table

| Signature | Dataset | Platform | Genes used | Samples | ROC-AUC | Interpretation |
|---|---|---:|---:|---:|---:|---|
| Ranked signature | Internal CV | Microarray | ~7 | 109 | ~0.95 | Best overall internal result |
| Ranked signature | GSE135485 | RNA-seq | 6 | 58 | 0.66 | Partial cross-platform transfer |
| Ranked signature | GSE25628 | Microarray | 2 | 14 | 0.77 | Good cross-study transfer |
| L1 sparse signature | Internal CV | Microarray | 3 | 109 | 0.91 | Strong sparse core |
| L1 sparse signature | GSE135485 | RNA-seq | 3 | 58 | 0.42 | Too compressed for cross-platform transfer |
| L1 sparse signature | GSE25628 | Microarray | 1 | 14 | 0.90 | Strong but fragile single-gene transfer |

---

## Signature Size Analysis

![Signature Size](results/figures/signature_size_curve.png)

Performance improves rapidly as genes are added, then reaches a plateau at approximately 7 genes.

### Interpretation

This indicates that the predictive signal is low-dimensional. Even though the original dataset contains tens of thousands of probes, most of the useful signal is concentrated in a small subset.

This is one of the key findings of the project:
- only a few genes are needed to capture most of the predictive information
- adding more genes eventually provides no meaningful gain

---

## Sparse L1 Signature

L1-regularized logistic regression was used to test whether the signal could be compressed further.

### Best sparse configuration
- ROC-AUC: approximately 0.91
- Non-zero genes: 3
- Selected genes:
  - ZNF24
  - HMGN3-AS1
  - ZNF568

### Interpretation

This confirms that the signal is genuinely sparse.

However, external validation shows that the 3-gene signature is less robust than the larger ~7-gene signature, especially across platforms. This suggests that part of the larger signature acts as stabilizing redundancy.

In other words:
- the 3-gene signature is more interpretable
- the 7-gene signature is more robust

---

## Signature Size Analysis

![Signature Size](results/figures/signature_size_curve.png)

To evaluate the number of genes required to capture the predictive signal, model performance was analyzed as a function of signature size.

### Results

- Performance improves as more genes are added
- A plateau is reached at approximately **7 genes**
- Adding additional genes does not improve ROC-AUC

### Interpretation

This indicates that the endometriosis signal is **low-dimensional** and can be captured by a small subset of genes.

Even though the original dataset contains tens of thousands of features, most of the predictive information is concentrated in a compact signature.

This result supports the use of minimal gene panels for downstream applications such as biomarker development.  

---

## External Validation

### GSE135485 (RNA-seq)
- Ranked signature ROC-AUC: 0.66
- Sparse L1 signature ROC-AUC: 0.42

This dataset is highly imbalanced, with 54 positives and only 4 negatives, so results should be interpreted carefully. Still, the comparison is informative: the larger signature transfers partially, while the minimal 3-gene signature fails.

### GSE25628 (microarray)
- Ranked signature ROC-AUC: 0.77
- Sparse L1 signature ROC-AUC: 0.90

For this dataset, only a subset of the signature was available after mapping across platforms. Even so, transfer was strong, showing that some genes, particularly ZNF24, carry a highly concentrated signal in this cohort.

### Overall interpretation

The project supports three conclusions:
1. Endometriosis-related signal can be captured by a compact signature.
2. The signal generalizes better across studies than across technologies.
3. Very aggressive compression improves interpretability but can reduce robustness.

---

## ROC Curves

![ROC Curves](results/figures/roc_curve_all.png)

The combined ROC plot shows:
- strong internal performance
- reduced but non-random transfer to RNA-seq
- stronger transfer to an independent microarray study

---

## Core Genes

The main interpretable signature includes genes such as:
- CTU2
- ZNF24
- NT5DC3
- HMGN3-AS1
- ZNF568
- C11orf54

The sparse L1 core is:
- ZNF24
- HMGN3-AS1
- ZNF568

These results suggest a distributed but compressible transcriptional signature rather than a single dominant biomarker.

---

## Project Structure

```text
.
├── data/
│   ├── raw/
│   └── processed/
├── scripts/
│   ├── build_dataset.py
│   ├── run_signature.py
│   ├── run_signature_l1.py
│   ├── validate_gse135485.py
│   ├── validate_gse25628.py
│   ├── plot_roc_all.py
│   └── plot_signature_size.py
├── src/
│   └── endometriosis_signature/
├── results/
│   └── figures/
└── README.md
```

---

## Usage

Install the package:

```bash
pip install -e .
```

Run the main pipeline:

```bash
python scripts/build_dataset.py
python scripts/run_signature.py
python scripts/run_signature_l1.py
python scripts/validate_gse135485.py
python scripts/validate_gse25628.py
python scripts/plot_signature_size.py
python scripts/plot_roc_all.py
```

---

## Data Availability

Raw data is not included in the repository due to size constraints.

Download from GEO:
- GSE51981
- GSE135485
- GSE25628

Place files in:

```text
data/raw/
```

Platform annotation files are also required for probe-to-gene mapping:
- GPL570
- GPL571

---

## Limitations

- External validation cohorts are small
- The RNA-seq validation set is highly imbalanced
- Cross-platform transfer remains difficult
- Biological interpretation is still preliminary
- No pathway enrichment analysis is included yet

---

## Future Work

Possible extensions include:
- pathway enrichment analysis
- larger and more balanced external validation cohorts
- nested cross-validation
- comparison with alternative linear baselines
- literature-guided biological interpretation of the sparse core genes

---

## Conclusion

This project shows that endometriosis-related transcriptional signal can be captured by a compact and stable gene signature.

The main takeaways are:
- the signal is low-dimensional
- a small signature performs nearly as well as the full model
- a 3-gene sparse core exists
- larger compact signatures are more robust for transfer
- generalization is stronger across studies than across platforms

Overall, the project illustrates both the promise and the limitations of machine learning for transcriptomic biomarker discovery.
