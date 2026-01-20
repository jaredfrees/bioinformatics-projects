library(hdf5r)
library(Seurat)
#library(SeuratData)
library(dplyr)
library(patchwork)
library(ggplot2)
library(hdf5r)
library(sctransform) ## sctransform method Replace need to run NormalizeData FindVarialbeFeature ScaleData with ScTransform

working_dir <- getwd()
setwd(working_dir)

raw_data_file_name <- 'S3d1_vpM_filtered_feature_bc_matrix.h5'
print(paste(working_dir, raw_data_file_name, sep=''))
d1vpm.data <- Read10X_h5(paste(working_dir, raw_data_file_name, sep=''))
d1vpm <- CreateSeuratObject(counts = d1vpm.data, project = "S3d1_vpM", min.cells = 3, min.features=200)


d1vpm[["percent.mt"]] <- PercentageFeatureSet(d1vpm, pattern = "^MT-")
VlnPlot(d1vpm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

rds_file_name <- "d1vpm.rds"
saveRDS(d1vpm, file = paste(working_dir, rds_file_name, sep=''))

getwd()
# Save last plot displayed
ggsave('d1_vpm_percent_mt.png', width=7, height=6, dpi = 600)

plot1 <- FeatureScatter(d1vpm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(d1vpm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(d1vpm, subset = nFeature_RNA > 10 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

head(d1vpm)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


