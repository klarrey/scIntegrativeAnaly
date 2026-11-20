
# load required libraries
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(harmony)
library(reticulate)

# install SeuratData datasets
InstallData("pbmcref")
InstallData("pbmcsca")


# load in the pbmc systematic comparative analysis dataset
obj <- LoadData("pbmcsca")
View(obj)
head(obj@meta.data)
dim(obj@meta.data)
# extract cells which meet the nFeature_RNA > 1000 condition
obj  <- subset(obj, nFeature_RNA > 1000)
# perform cell annotations
obj <- RunAzimuth(obj, reference = "pbmcref")
# currently, the object has two layers in the RNA assay: counts, and data
obj

# Split the layers
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
obj


obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# visualize the results of the standard analysis without integration
obj <- FindNeighbors(obj, dims=1:30, reduction = 'pca')
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = 'pca', reduction.name = 'umap.unintegrated')
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(obj, reduction = 'umap.unintegrated', group.by=c('Method','predicted.celltype.l2'))

# Perform streamlined (one-line) integrative analysis
obj_cca <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.cca',
  verbose = FALSE)

obj_rpca <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.rpca',
  verbose = FALSE)

obj_harmony <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = 'harmony',
  verbose = FALSE)

obj_fastmnn <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = 'integrated.mnn',
  verbose = FALSE)

obj_scvi <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  new.reduction = 'integrated.scvi',
  conda_env = '../miniconda3/envs/scvi-env', verbose = FALSE)

# clustering for ccaintegration
obj_cca <- FindNeighbors(obj_cca, reduction = 'integrated.cca', dims = 1:30)
obj_cca <- FindClusters(obj_cca,resolution = 2, cluster.name = 'cca_clusters')

# visualization for ccaintegration
obj_cca <- RunUMAP(obj_cca, reduction = "integrated.cca", dims = 1:30, reduction.name = 'umap.cca')
p1 <- DimPlot(
  obj_cca, reduction = "umap.cca",
  group.by = c("Method", "predicted.celltype.l2", "cca_clusters"),
  combine = FALSE, label.size = 2) 

obj_scvi <- FindNeighbors(obj_scvi, reduction = 'integrated.scvi', dims = 1:30)
obj_scvi <- FindClusters(obj,resolution = 2, cluster.name = 'scvi_clusters')

obj_scvi <- RunUMAP(obj_scvi, reduction = "integrated.scvi", dims = 1:30, reduction.name = 'umap.scvi')
p2 <- DimPlot(
  obj_scvi, reduction = "umap.scvi",
  group.by = c("Method", "predicted.celltype.l2", "scvi_clusters"),
  combine = FALSE, label.size = 2)

wrap_plots(c(p1, p2), ncol = 2, byrow = F)