# Pscore
Pscore is a comprehensive toolkit for survival analysis in prostate cancer research, integrating deep learning models with clinical and genomic data to predict patient outcomes.

## Overview
This project combines multiple data modalities (WSI images, gene expression, clinical features) to create an integrated prognostic score for cancer patients. The core components include:

- Deep learning survival model (DeepSurv)
- Differential gene expression analysis
- Gene Set Enrichment Analysis (GSEA)
- Multi-modal data integration and training
## Components
### DeepSurv
A deep learning framework for survival analysis with the following features:

- Neural network architecture for risk prediction
- Negative log-likelihood loss function
- Automatic weighted loss for multi-task learning
- Cross-validation and ensemble prediction
### Scoring System
The project implements multiple scoring approaches:

- Gene expression score (Gscore)
- Clinical feature score (Cscore)
- WSI-based score (Wscore)
- Integrated prognostic score (Pscore)
### Analysis Tools
- DEGs.R : Differential gene expression analysis
- GSEA.r : Gene Set Enrichment Analysis
- TrainScore.r : Training and evaluation of scoring models
- TestScore.R : Testing and validation of scoring models
- compare.R : Compare Pscore with PMID: 39296545(https://pubmed.ncbi.nlm.nih.gov/39296545/) and PMID: 36439479(https://pubmed.ncbi.nlm.nih.gov/36439479/)

## Requirements
### Python Dependencies
- PyTorch
- NumPy
- Pandas
- h5py
### R Dependencies
- tidyverse
- survival
- limma
- edgeR
- ggplot2
- survminer
## Acknowledgements
Special thanks to the following contributors who made this project possible:

- czifan (czifan@pku.edu.cn, https://github.com/czifan/DeepSurv.pytorch) for the excellent DeepSurv implementation
- Bowen Zheng for developing the AutomaticWeightedLoss component
