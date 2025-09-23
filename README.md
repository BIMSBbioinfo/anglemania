<img src="man/figures/graphical_abstract.png" align="center" alt="logo" width="2000" style = "border: none; float: center ;">

# anglemania
## Overview
**anglemania** is a new approach to the integration of scRNA-seq (and, potentially, others sc-omics) from **similar** biological entities.
The novelty, as well as the cornerstone, of the proposed approach, is to use the conservation of [angles](https://arxiv.org/abs/1306.0256) between gene pairs across an assembly of datasets to be integrated.
anglemania extracts genes from gene pairs that exhibit invariant and biologically meaningful relationships across different experiments. Those genes can subsequently be used as the basis for integration algorithms such as [Seurat CCA integration](https://satijalab.org/seurat/reference/ccaintegration), [SCVI](https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html) and others.
Checkout our [website](https://bioinformatics.mdc-berlin.de/anglemania).

|  |  | 
| - | - |
| Github | [![Github](https://github.com/BIMSBbioinfo/anglemania/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BIMSBbioinfo/anglemania/actions/workflows/R-CMD-check.yaml) |
| Bioc Release | [![Bioc Release](https://bioconductor.org/shields/years-in-bioc/anglemania.svg)](https://bioconductor.org/packages/anglemania)




## Installation

You can install **anglemania** from Bioconductor using the following commands:
**once available on Bioconductor:**

```r
# Install BiocManager if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install anglemania
BiocManager::install("anglemania")
```

For the development version:

```r
BiocManager::install("anglemania", version="devel")
```


## Documentation

Comprehensive documentation is available on our [anglemania website](https://bioinformatics.mdc-berlin.de/anglemania)
Or visit the [Bioconductor package page](https://bioconductor.org/packages/anglemania).

## Getting Help

If you encounter issues or have questions:

- Submit issues on [GitHub Issues](https://github.com/BIMSBbioinfo/anglemania/issues).

## License

Artistic License/GPL

## Citation

If you use **anglemania** in your research, please cite:

> Aaron Kollotzek, Vedran Franke, Artem Baranovskii, Altuna Akalin(2024). *anglemania: Feature Extraction for scRNA-seq Dataset Integration. R package version 0.99.6




