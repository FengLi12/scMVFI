A framework for scRNA-seq data clustering based on multi-view feature integration
##Overview

scMVFI is built upon a single-cell feature extraction framework that utilizes autoencoder networks and multi-view feature fusion. In this framework, the autoencoder network is employed to reconstruct gene expression values from scRNA-Seq data, alleviating the issue of dropout. Additionally, intrinsic entropy is utilized for feature pre-selection, and multiple similarity estimates from various feature sets are integrated from different perspectives to construct a comprehensive cell-to-cell similarity network.


## Requirement:

- `python = 3.9.12`
- `numpy = 1.15.4`
- `pandas = 0.23.4`
- `tensorflow-cpu = 2.9.1`
- `R = 4.3.1`
- `umap=0.2.10.0`
- `parallelDist=0.2.6`
- `irlba=2.3.5.1`
- `loe=1.1`




