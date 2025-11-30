# ATRT-Subgroup-Classifier 
## Reproducing analysis + Machine learning model

### Project overview

Atypical teratoid/rhabdoid tumors (ATRTs) are aggressive pediatric brain tumors.
They represent the most common malignant brain tumors that manifests in infancy.
ATRTs can be divided into three main subgroups (TYR, SHH and MYC) that differs in their gene expression profile, epigenetic signature and clinical behaviour.

This project replicates the unsupervised hierarchical cluster analysis of gene expression profiles
originally published by Johann, Pascal D. et al. 2016 to distinguish the three molecular subgroups of ATRTs.

Using RNA microarray data, I developed a machine learning pipeline to:

1. Reconstruct the 3 distinct subgroups using unsupervised clustering.

2. Train an ensemble of classifiers (RF, XGBoost, GLMNET, PLS-DA).

3. Classified two external independent datasets.

### Project workflow
This analysis is fully automated using a Snakemake pipeline to ensure reproducibility from start to finish.
![Snakemake-dags](/results/figures/dag.png)

### Quick start - Installation & Execution

#### Option A
```bash
# 1. Clone
git clone https://github.com/vini-rios/atrt
cd atrt

# 2. Setup environment
bash scripts/setup_project.sh

# 3. Run analysis
snakemake --cores all --use-conda
```
#### Option B - Docker
```bash
# 1. Clone
git clone https://github.com/vini-rios/atrt
cd atrt

# 2. Setup environment
bash scripts/setup_project.sh

# 3. Run analysis
snakemake --cores all --use-conda
```
### Requirements

* Linux / macOS (or Docker)

* R >= 4.0

* Snakemake >= 6.0

* Conda (for environment.yml) or Docker

* Python 3.8+ for some pipeline steps (Snakemake)