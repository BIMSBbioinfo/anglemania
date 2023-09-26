###############################################################
### AUTHOR:ab
### DATE: 10.08.2023
### INFO: Script to load all necessary libraries and functions 
###############################################################

## libraries ---- temporary
library(tidyverse)
library(tidygraph)
library(ggraph)
library(janitor)
library(tictoc)
library(data.table)
library(Matrix)
library(Seurat)
library(parallelDist)
library(digest)
##
library(progressr)
handlers(global = TRUE)
handlers("progress")
##
library(Rcpp)
sourceCpp("./src/dist2mat.cpp")

##
setDTthreads(16)
## functions
source("./src/anglemanise.R")
source("./src/factorise.R")
source("./src/flippitise.R")
source("./src/F_main.R")
source("./src/F_auxiliary.R")
source("./src/F_stats.R")