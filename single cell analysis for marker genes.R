#=======================================================================================#
# The usages of all packages:                                                           #
#      dplyr: provide a flexible grammar of data manipulation                           #
#     export: export the R plot into PPT, Word, jpg, png, tif, pdf, etc. table to excel #
#     Seurat: create an object for single-cell expression data to reduce the dimensions #
#    cowplot: supply a simple code for easy theme change while plotting                 #
#  tidyverse: assemble the packages of data manipulation and visualization              #
# reticulate: supply the R interface to Python for UMAP clustering                      #
#=======================================================================================#


##################
# Load libraries #
##################
# Clear and garbage
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
# Load packages
{
  library(openxlsx)
  library(Seurat)
  library(dplyr)
  library(export)
  library(cowplot)
  library(tidyverse)
  library(reticulate)
}

########################################
# Load Data and Create a Seurat Object #
########################################
# Load Seurat Object and Data
if(!file.exists("sce.Rdata")){
  counts <- read.xlsx("Table S1. Single cell gene expression profile.xlsx", rowNames = TRUE)
  sce <- CreateSeuratObject(
    counts = as.matrix(counts),
    project = "sc-RNA seq",
    min.cells = 12, 
    min.features = 200,
    assay="RNA")
  save(counts, sce, file = "sce.Rdata")
}
load("sce.Rdata")
dim(sce)
# [1] 13877   239

######################
# Data Normalization #
######################
# Normalize the data, log-transform the result and mutiply this by a scale factor
sce <- NormalizeData(
  object = sce,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

#########################
# Detect Variable Genes #
#########################
# Calculate highly variable genes and focuse on these for downstream analysis: Calculate the average expression 
# and dispersion for each gene, place these genes into bins and calculate a z-score for dispersion to control for 
# the relationship between variablity and average expression
sce <- FindVariableFeatures(
  object = sce,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  x.low.cutoff = 0.0125,
  x.high.cutoff = 3,
  y.cutoff = 0.5)
head(x = HVFInfo(object = sce))

#############################
# Data Scaling and Removing #
#############################
# Scale and center the data to regress out the noise-like siganal and nCounts_scRNA
sce <- ScaleData(object = sce, vars.to.regress = "nCount_RNA")

########################################
# Perform linear dimensional reduction #
########################################
# Perform the PCA reduction for scaled data
sce <- RunPCA(object = sce, verbose = FALSE)
# Examine and visualize the PCA result
DimPlot(object = sce,reduction = "pca")
# Scatter plot across single cells without geneplot
FeatureScatter(object = sce, feature1 = "Jag2", feature2 = "PC_1")
FeatureScatter(object = sce, feature1 = "Jag2", feature2 = "Jag1")
# Scatter plot across individual features without cellplot
CellScatter(
  object = sce,
  cell1 = "VB001",
  cell2 = "VB002")
# Divide the counts into non-variable and variable count
VariableFeaturePlot(object = sce)
# Violin nad Ridge plots
VlnPlot(object = sce,features = c("Notch1", "Notch2", "Notch3"))
# Export the graph to PPT
# graph2ppt(file = "Violin plot.pptx", width=7, height=5)
RidgePlot(object = sce,feature = c("Notch1", "Notch2", "Notch3"))
# Heatmaps
DimHeatmap(object = sce, reduction = "pca", cells = 200, balanced = TRUE)

###################
# Cell Clustering #
###################
# Calculate the k-nearest neighbors and construct the SNN graph
sce <- FindNeighbors(sce, reduction = "pca")
sce <- FindClusters(sce, resolution = 0.5, algorithm = 1)

############################################
# Perform non-linear dimensional reduction #
############################################
# Visualize and explore the datasets using tSNE
sce <- RunTSNE(object = sce)
DimPlot(object = sce, reduction = "tsne", do.label = TRUE, pt.size = 2, order = 7)

#######################################
# Find Differentially Expressed genes #
#######################################
# Decide markers for every cluster compared to all other cells
VlnPlot(object = sce, features = c("Notch1", "Notch2", "Notch3", "Jag2", "Dll1"))
FeaturePlot(object = sce, feature = c("Cd93", "Flt1", "Igfbp7"), 
            cols=c("grey", "red"), reduction="tsne", pt.size = 2)
FeaturePlot(object = sce, feature = c("Dll1", "Jag2", "Notch1", "Notch2", "Notch3"),
            cols=c("grey", "red"), reduction="tsne", pt.size = 2)

# Save Data
save(sce, file = "sce_analysed.Rdata")

#============================#
#       Musician: Resonance  #
#           Date: 2019/10/14 #
# Revised author: Resonance  #
#           Time: 2019/10/14 #
#============================#