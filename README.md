<img src="800px-Tissot_The_Gathering_of_the_Manna_(color).jpg" align="left" alt="logo" width="300" style = "border: none; float: center ;">

# anglemana
The repository carries the development version of the conceived "anglemana" R package alongside the ongoing tasks and obstacles. 

## Introduction
The goal of this project is to develop a new approach to the integration of scRNA-seq (and, potentially, others sc-omics?) from **similar** biological entities.
The novelty, as well as the cornerstone, of the proposed approach, is to use the conservation of [angles](https://arxiv.org/abs/1306.0256) between gene pairs across an assembly of datasets to be integrated. 
The approach relies on preprocessing of the scRNA-seq data with [SEACELL](https://www.nature.com/articles/s41587-023-01716-9), which is used as an input.
Currently, the project works as a methodology for feature selection, while integration is done by CCA + mNN as implemented in Seurat. 

## Status
We have achieved stable performance across multiple neuroblastoma datasets, i.e. conserved topology was recovered.
A script ```anglemanise.R``` contains a master function that runs extracts the conserved angles between datasets. 
A template run can be found in ```playground/templateRun.R```

## TODOs
- Implement multimodality testing during the selection of critical angles (```src/F_stats.R```).
- Write an algorithm that iteratively builds a graph between nodes using conserved angles as edges and stops when a certain connectivity criterion is reached, e.g. [centrality_communicability()](https://tidygraph.data-imaginist.com/reference/centrality.html)
- Conceive an approach to account for a vector flip (flippity), i.e. a situation where the angle remains the same while relative positions within a gene pair are invested. 
  - PCA reference is insufficient.
