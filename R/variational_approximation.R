suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(data.table)
    library(Seurat)
    library(countsplit)
    library(Matrix)
})


split_matrix = function(X, epsilon = 0.5){

    Xtrain <- apply(X, 2, function(u) rbinom(n = length(u),
    size = as.integer(u), prob = epsilon))
    rownames(Xtrain) <- rownames(X)
    colnames(Xtrain) <- colnames(X)
    Xtest <- X - Xtrain
    rownames(Xtest) <- rownames(X)
    colnames(Xtest) <- colnames(X)
    return(list(train = Xtrain, test = Xtest))
}






# ------------------------------------------------------------------------------------- #
#' Find consensus variable features from a list of Seurat objects
#'
#' This function takes an seu object and processes all samples. It normalizes the data, scales it, 
#' finds variable features, and then produces a data frame of all variable features that are present 
#' in more than 4 samples.
#'
#' @param seu_list a list of seu objects
#' @param nfeatures the number of variable features to find in each sample
#' @param loess_span the span of the loess function for finding variable features
#' @param n_samp minimal number of samples in which the feature must be present
#'
#' @return a data frame of features that are present in more than 4 samples
find_consensus_variable_features = function(
  seu_list, 
  nfeatures = 1000, 
  loess_span = 0.4,
  n_samp = 5,
  filter_dataset = FALSE
  
) {
   
  # Extract the sample names
  if(is.null(names(seu_list)))
    stop("seu_list has to be a named list")
  sample_names = names(seu_list)
  seu_list = lapply(seu_list, function(x){
    DefaultAssay(x) = "RNA"
    x
  })
  # For each sample, normalize, scale, and find variable features
  message("Variable features ...")
  feature_list = lapply(setNames(sample_names, sample_names), function(x){
    message(x)
    seu_list[[x]] %>%
      NormalizeData(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      FindVariableFeatures(verbose = FALSE, nfeatures = nfeatures, loess.span = loess_span) %>%
      VariableFeatures()
  })

  # Create a data frame of the variable features for each sample
  message("Consensus features ...")
  dfeature = lapply(names(feature_list), function(x){
    data.frame(sample_name = x, features = feature_list[[x]])
  }) %>%
    bind_rows() %>%
    left_join(lapply(seu_list, function(x){
        x@meta.data %>%
          dplyr::select(sample_name, dataset) %>%
          distinct() 
      }) %>% bind_rows()
    )

    ### Checks whether the feature is present in multiple datasets
    if(filter_dataset){
      message("using dataset")
      filtered_feature = dfeature %>% 
        group_by(features) %>%
        mutate(n_dataset = length(unique(dataset))) %>%
        ungroup() %>%
        mutate(n_datasets = length(unique(dataset))) %>%
        filter(n_datasets == n_dataset) %>%
        ungroup() %>%
        group_by(features) %>%
        mutate(n = n()) %>%
        ungroup() %>%
        filter(n >= n_samp) %>%
        filter(!duplicated(features))
    }else{
      message("using sample")
      filtered_feature = dfeature %>% 
        group_by(features) %>%
        mutate(n = n()) %>%
        filter(n >= n_samp) %>%
        ungroup() %>%
        filter(!duplicated(features))
    }
    
  
  return(filtered_feature)
}

# ------------------------------------------------------------------------------------- #
#' Integrate a list of Seurat objects
#'
#' This function takes a list of Seurat objects and then normalizes and scales the data.
#' It uses prespecified features, or finds the necessary features
#' for visualization and clustering.
#'
#' @param seu_list  a list of Seurat objects objects
#' @param dfeature data frame of features to use for anchor identification
#' @param reduction the dimensionality reduction method to use
#' @param nfeatures the number of features to use for anchor identification
#' @param dims the dimensions to consider in various operations
#' @param k_filter, k_anchor, k_score parameters for FindIntegrationAnchors function
#' @param k_weight parameter for IntegrateData function
#' @param resolution resolution for FindClusters function
#'
#' @return a Seurat object after all processing and integration
Integrate_Seurat_List_V5 = cacheFile(path_RDS) %@% function(
  seu_list, 
  features = NULL, 
  reduction = "cca", 
  dims = 1:10,
  k_filter = 10, 
  k_anchor = 10, 
  k_score = 10, 
  k_weight = 10, 
  resolution = seq(0.4, 2, 0.2),
  integration_features = 2000,
  filter_cell_number = 3,
  filter_cnts_value  = 0
){
  options(future.globals.maxSize = 3 * 1024^3)
  
  # Normalize and scale the data for each sample
  sample_names = names(seu_list)
  seu_list = lapply(setNames(sample_names, sample_names), function(x){
    message(x)
    seu = seu_list[[x]] %>%
      NormalizeData(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(verbose = FALSE, features = features, npc = 15)
    LayerData(seu, "scale.data") = NULL
    seu
  })
    
  if(is.null(features)){
    message("Features not supplied, finding anchor features ...")
    features = SelectIntegrationFeatures(
      object.list = seu_list, 
      nfeatures = integration_features, 
      verbose=FALSE
    )
  }

  # Find integration anchors and integrate the data
  message("Merge ...")
  
  seu_merged =  merge(seu_list[[1]], seu_list[-1], merge.dr = TRUE)
  seu_merged = ScaleData(seu_merged)

  message("Integrate layers ...")
  seu_combined = IntegrateLayers(
    seu_merged, 
    method   = "CCAIntegration", 
    layers   = Layers(seu_merged, "data"), 
    assay    = "RNA",
    features = features,
    dims     = dims,
    k.filter  = k_filter, 
    k.anchor  = k_anchor, 
    k.score   = k_score,
    k.weight  = k_weight
  )
  seu_combined = JoinLayers(seu_combined)
 
  # Run the standard workflow for visualization and clustering
  message("Processing ...")
  seu_combined = ScaleData(seu_combined, verbose = FALSE) 
  seu_combined = RunUMAP(
      seu_combined,
      reduction = "integrated.dr",
      dims = 1:ncol(Embeddings(seu_combined$integrated.dr)),
      return.model = TRUE
  )
  set.seed(12345)
  seu_combined = FindNeighbors(seu_combined, reduction = "umap", dims = 1:2, verbose = TRUE)
  seu_combined = FindClusters(seu_combined, resolution = resolution)
  DefaultAssay(seu_combined) = "RNA"
  
  return(seu_combined)
}


# ------------------------------------------------------------------------------------- #
#' Integrate a list of Seurat objects
#'
#' This function takes a list of Seurat objects and then normalizes and scales the data.
#' It uses prespecified features, or finds the necessary features
#' for visualization and clustering.
#'
#' @param seu_list  a list of Seurat objects objects
#' @param dfeature data frame of features to use for anchor identification
#' @param reduction the dimensionality reduction method to use
#' @param nfeatures the number of features to use for anchor identification
#' @param dims the dimensions to consider in various operations
#' @param k_filter, k_anchor, k_score parameters for FindIntegrationAnchors function
#' @param k_weight parameter for IntegrateData function
#' @param resolution resolution for FindClusters function
#'
#' @return a Seurat object after all processing and integration
Integrate_Seurat_List_V4 = cacheFile(path_RDS) %@% function(
  seu_list, 
  features = NULL, 
  reduction = "cca", 
  dims = 1:10,
  k_filter = 10, 
  k_anchor = 10, 
  k_score = 10, 
  k_weight = 10, 
  resolution = seq(0.4, 2, 0.2),
  integration_features = 2000,
  filter_cell_number = 3,
  filter_cnts_value  = 0
){
  options(future.globals.maxSize = 3 * 1024^3)
  
  # Normalize and scale the data for each sample
  sample_names = names(seu_list)
  seu_list = lapply(setNames(sample_names, sample_names), function(x){
    message(x)
    seu = seu_list[[x]] %>%
      NormalizeData(verbose = FALSE) %>%
      ScaleData(verbose = FALSE)
  })
    
  if(is.null(features)){
    message("Features not supplied, finding anchor features ...")
    features = SelectIntegrationFeatures(
      object.list = seu_list, 
      nfeatures = integration_features, 
      verbose=FALSE
    )
  }

  # Find integration anchors and integrate the data
  message("Find Anchors ...")
  anchors = FindIntegrationAnchors(
      object.list = seu_list, 
      anchor.features = features, 
      verbose = FALSE,
      dims = dims,
      k.filter  = k_filter, 
      k.anchor  = k_anchor, 
      k.score   = k_score,
      reduction = reduction
  )
  features_intersect = Reduce(function(x,y)intersect(x,y), lapply(seu_list, rownames))
  seu_combined = IntegrateData(
    anchorset = anchors, 
    k.weight = k_weight, 
    features.to.integrate = features_intersect
  )
  DefaultAssay(seu_combined) = "integrated"
  
  # Run the standard workflow for visualization and clustering
  message("Processing ...")
  seu_combined = ScaleData(seu_combined, verbose = FALSE)
  seu_combined = RunPCA(seu_combined, npcs = 30, verbose = FALSE)
  seu_combined = RunUMAP(
      seu_combined,
      reduction = "pca",
      dims = 1:30,
      umap.method = "uwot-learn",
      return.model = TRUE
  )
  set.seed(12345)
  seu_combined = FindNeighbors(seu_combined, reduction = "umap", dims = 1:2, verbose = TRUE)
  seu_combined = FindClusters(seu_combined, resolution = resolution)
  DefaultAssay(seu_combined) = "RNA"
  
  return(seu_combined)
}

# ------------------------------------------------------------ #

FindMarkers_All_Clusters = cacheFile(path_RDS) %@% function(
  seu,
  cluster,
  min.pct = 0.1,
  assay = "integrated",
  only.pos = TRUE,
  logfc.threshold = 0.5,
  test_use = "wilcox"
){
  Idents(seu) = seu@meta.data[[cluster]]
  FindAllMarkers(
    seu, 
    only.pos        = only.pos, 
    logfc.threshold = logfc.threshold, 
    assay           = assay, 
    min.pct         = min.pct,
    test.use        = test_use
  )
}

# ---------------------------------------------------------------- #
#' Optimize Integration Parameters for Single Cell Data Integration
#'
#' This function optimizes integration parameters for single cell data by
#' iterating over a grid of parameter values, performing data integration,
#' and saving the results. It is designed to work with Seurat objects.
#'
#' @param seu A Seurat object containing single cell data to be integrated.
#' @param outpath_param_optimize The directory path where the optimization
#' results will be saved. The function will create this directory if it does
#' not exist.
#' @param integration_params A data frame of integration parameters to iterate
#' over, specifically the number of features and the number of samples. Default
#' values range for features (nfeatures) are 500, 1000, 2000, and 4000, and for
#' number of samples (nsamp) are 7 to 13.
#' @param ncores The number of cores to use for parallel processing. Default
#' is 16.
#' @param split_by The column name in the metadata of the Seurat object that
#' indicates how to split the object for integration. Default is "sample_name".
#'
#' @return A list of Seurat objects resulting from the integration process with
#' different parameters.
#' @import doMC
#' @import foreach
#' @examples
#' # Example usage (assuming 'seu' is a Seurat object):
#' results = Optimize_Integration_Params(
#'   seu,
#'   "path/to/output",
#'   ncores = 4,
#'   split_by = "sample_group"
#' )
#' @export
#'
#' ---------------------------------------------------------------- #
optimize_integration_params = function(
  seu,
  outpath_param_optimize,
  integration_params = expand.grid(
    nfeatures = c(500, 1000, 2000, 4000),
    nsamp     = c(7:13)
  ),
  ncores = 16,
  split_by = "sample_name"
){
  # Create the output directory if it doesn't exist
  dir.create(outpath_param_optimize, showWarnings=FALSE)

  # Load parallel processing library and register cores
  library(doMC)
  library(patchwork)
  registerDoMC(ncores)
  
  # Source the script for single cell data integration
  source.pro("Integrate_Single_Cell.R")

  # Iterate over each combination of integration parameters
  lres = foreach(i = 1:nrow(integration_params), .errorhandling="pass")%dopar%{
    message(i) # Output the iteration index

    source.pro("Integrate_Single_Cell.R") # Re-source the integration script (necessary if it defines functions)
    
    # Extract the current parameters
    nsamp     = integration_params$nsamp[i]
    nfeatures = integration_params$nfeatures[i]
    
    # Set the default assay
    DefaultAssay(seu) = "RNA"
    
    # Split the Seurat object by the specified metadata column
    seu_train_list = SplitObject(seu, split.by = split_by)
    
    # Skip if the number of samples is greater than the number of split objects
    if(nsamp > length(seu_train_list))
      next()

    # Find consensus variable features based on the current parameters
    consensus_features = find_consensus_variable_features(seu_train_list, nfeatures=nfeatures, n_samp = nsamp)
    
    # Skip if the number of consensus features is less than 100
    if(nrow(consensus_features) < 100)
      next()
      
    # Integrate the Seurat objects based on the consensus features
    seu_combined = Integrate_Seurat_List_V4(
      seu_train_list, 
      features = consensus_features$features
    )

    # Generate UMAP plots for the integrated object, grouping by sample name and dataset
    p1 = DimPlot(seu_combined, reduction = "umap", group.by = "sample_name")
    p2 = DimPlot(seu_combined, reduction = "umap", group.by = "dataset")
  
    # Combine and save the plots
    gg = p1 + p2 + plot_layout(nrow=1)
    ggsave(file.path(outpath_param_optimize, paste0("umap-combined-SampleName_Cluster_Annot","_nfeat",nfeatures,"_nsamp",nsamp,".png")), gg, width = 15, height = 5)
    
    # Construct a name for the parameter combination
    params_name = paste("nfeat",nfeatures,"nsamp",nsamp, sep="_")
    
    # Return the integrated Seurat object
    return(seu_combined)
  }
  invisible(return(lres))
}
