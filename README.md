<img src="800px-Tissot_The_Gathering_of_the_Manna_(color).jpg" align="left" alt="logo" width="300" style = "border: none; float: center ;">

# anglemania
The repository carries the development version of the conceived "anglemania" R package alongside the ongoing tasks and obstacles. 

## Introduction
The goal of this project is to develop a new approach to the integration of scRNA-seq (and, potentially, others sc-omics?) from **similar** biological entities.
The novelty, as well as the cornerstone, of the proposed approach, is to use the conservation of [angles](https://arxiv.org/abs/1306.0256) between gene pairs across an assembly of datasets to be integrated. 
The approach relies on preprocessing the scRNA-seq data with [SEACell](https://www.nature.com/articles/s41587-023-01716-9), which is used as an input.
Currently, the project works as a methodology for feature selection, while integration is done by CCA as implemented in Seurat. 

## Status
We have achieved stable performance in Jansky and Dong neuroblastoma datasets. However, the manual tuning is required for both.
Other datasets need to be tested. An improvement over the *status quo* is expected. A few ideas are: the adjustment of the integration order, a smart approach to feature selection from conserved angles, and scaling of the integration feature importance by angle conservation, strict expression filtering before extraction of gene-gene angles, e.g. some angles show a conserved relationship, although the expression is negligibly small, which invalidates their usefulness.

Currently, all used functions are documented and ready to be roxygenised. Technically, the package can be built and tested as soon as the latter is done.


## Ideas

Following [google doc](https://docs.google.com/document/d/10TEWmnfBOlW7SGFl70eb_26Z-kTeOK-bcNg6WgZJMsk/edit?usp=sharing)

<!-- badges: start -->
  [![R-CMD-check](https://github.com/BIMSBbioinfo/anglemania/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BIMSBbioinfo/anglemania/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
