# Set working directory
setwd("~/Downloads/TUMOR/patient1")
setwd("/Users/SHRITHIKA/Desktop/shrithika-csl/TUMOR/patient1")

# Load required packages
library(Seurat)
library(sceasy)
library(SeuratDisk)
library(remotes)
library(umap)
library(dplyr)
library(BiocManager)
library(BiocFileCache)

# Load data from Loom format
obj <- LoadLoom("1ad609b5-8c9f-4741-8a89-df03b6808b6d.seurat.1000x1000.loom")
print(obj)

# -------- 1. Quality Control (QC) -------- #
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plot to check data quality
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1)

# -------- 2. Data Filtering -------- #
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# -------- 3. Normalize Data -------- #
obj <- NormalizeData(obj)

# -------- 4. Identify Highly Variable Genes -------- #
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(obj), 10)

# Plot Variable Features
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# -------- 5. Scale Data -------- #
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)

# -------- 6. Principal Component Analysis (PCA) -------- #
obj <- RunPCA(obj, features = VariableFeatures(obj))
VizDimLoadings(obj, dims = 1:2, reduction = "pca")
DimPlot(obj, reduction = "pca") + NoLegend()
ElbowPlot(obj)

# -------- 7. Clustering -------- #
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.1)

# Assign cluster identities
new.cluster.ids <- c(
  "Astrocyte", "Excitatory/Inhibitory/OPC/Tumor", "Microglial", "T cell",
  "Unassigned 1", "Oligodendrocyte", "Endothelial", "Unassigned 2"
)
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

# -------- 8. UMAP Visualization -------- #
obj <- RunUMAP(obj, dims = 1:15)

# Define colors for visualization
custom_colors <- c(
  "Astrocyte" = "#FFD580", "Excitatory/Inhibitory/OPC/Tumor" = "#CCCCFF",
  "Microglial" = "#C27BA0", "T cell" = "#90EE90", "Unassigned 1" = "#D3D3D3",
  "Oligodendrocyte" = "#93CCEA", "Endothelial"= "#DAA06D", "Unassigned 2" = "#A9A9A9"
)

# Save UMAP plot
tiff("/Users/SHRITHIKA/Downloads/anubha/Figure/Supplementry_fig/Sup_fig1/umap_p1.tiff",
     units="in", width=9, height=5, res=300)
DimPlot(obj, reduction = "umap") + 
  scale_color_manual(values = custom_colors) +
  theme_minimal()
dev.off()
