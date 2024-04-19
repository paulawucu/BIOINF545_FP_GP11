#####ADDING LIBRARIES####
# Install necessary packages if not already installed
packages_to_install <- c(
  "BiocManager", "Seurat", "Matrix", "ggplot2", "dplyr", "devtools", "R.utils", 
  "harmony", "remotes", "SeuratDisk","SeuratDisk"
)

# Check if each package is installed and install if not
for (package in packages_to_install) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}

# Install specific Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'lme4', 
  'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 
  'HDF5Array', 'terra', 'ggrastr', 'AnnotationDbi', 'org.Hs.eg.db', 'DDRTree', 
  'pheatmap', 'DESeq2'
))

# Install packages from GitHub
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github('cole-trapnell-lab/monocle3')
remotes::install_github("satijalab/Azimuth")

# Load libraries
library(Seurat)
library(harmony)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Matrix)
library(tidyverse)
library(monocle3)
library(SeuratWrappers)
library(gridExtra)
library(Azimuth)
library(limma)
library(DESeq2)
library(SeuratDisk)
library(SeuratData)

#####UPLOADING DATA#####

raw_data_ys <- Read10X(data.dir = "matrix_files_ys_meyloid_progen_cells")

metadata_ys <- read.csv("metadata_ys.csv")

rownames(metadata_ys) <- metadata_ys$X
  
#####MAKING SEURAT OBJECT#####

sobj_ys <- CreateSeuratObject(counts = raw_data_ys, meta.data = metadata_ys)

view(sobj_ys@meta.data)

sobj_ys$mitopercent <- PercentageFeatureSet(sobj_ys, pattern = "^MT-")

sobj_ys_filtered <- subset(sobj_ys, subset = nCount_RNA > 800 &
                             nFeature_RNA >500 &
                             mitopercent < 10)

#Subsetting


#### MARIANA ADDED TO SOBJ_YS TO GET PCA####
sobj_ys_filtered$orig.ident <- "yolk_sac"

# Normalize data
sobj_ys_filtered <- NormalizeData(sobj_ys_filtered)

# Find variable features
sobj_ys_filtered <- FindVariableFeatures(sobj_ys_filtered)

# Scale data
sobj_ys_filtered <- ScaleData(sobj_ys_filtered)

# Compute PCA
sobj_ys_filtered <- RunPCA(sobj_ys_filtered, verbose = FALSE)

#Finding Neighboors
sobj_ys_filtered <- FindNeighbors(sobj_ys_filtered,dims = 1:30)

#clustering
sobj_ys_filtered <- FindClusters(sobj_ys_filtered, resolution = .9)

# Now run UMAP
sobj_ys_filtered <- RunUMAP(sobj_ys_filtered, dims = 1:30, n.neighbors = 50)

view(sobj_ys_filtered@meta.data)

A1 <- DimPlot(sobj_ys_filtered, reduction = 'umap', group.by = 'fetal.ids', label = T)
  
A2 <- DimPlot(sobj_ys_filtered, reduction = 'umap', group.by = 'seurat_clusters', label = T)

A1|A2

ElbowPlot(sobj_ys_filtered)

#####MARIANA ADDED: INTEGRATION #####

##Integration Step using Harmony: Christian's Method
#sobj_ys_integrated <- IntegrateLayers(object = sobj_ys_filtered, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",)
#sobj_ys_integrated <- FindNeighbors(sobj_ys_integrated, reduction = "harmony", dims = 1:8)
#sobj_ys_integrated <- FindClusters(sobj_ys_integrated, resolution = 0.3, cluster.name = "harmony_cluster")
#sobj_ys_integrated <- RunUMAP(sobj_ys_integrated, reduction = "harmony", dims = 1:8, reduction.name = "harmony.umap")

#I'm just splitting by component here since they were run on different lanes and that seems to be the batch effect that's overtaking the data right now

obj.list <- SplitObject(sobj_ys_filtered, split.by = 'fetal.ids')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

dim(sobj_ys_filtered)

# Loop through obj.list and check dimensions
for (i in seq_along(obj.list)) {
  cat("Object", i, "dimensions:", dim(obj.list[[i]]), "\n")
}

#removing object 3 that's less than 30 cells
obj.list <- obj.list[-3]

# Loop through obj.list and check dimensions: Checking if object 3 is removed 
for (i in seq_along(obj.list)) {
  cat("Object", i, "dimensions:", dim(obj.list[[i]]), "\n")
}
# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

#finding anchors 
# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)
# integrate data
sobj_ys_filtered_integrated <- IntegrateData(anchorset = anchors,k.weight = 50)

#Increasing k.anchor will result in selecting more anchor features from each dataset, potentially improving the integration accuracy, but also increasing computational cost.
#Lowering k.weight reduces the number of anchor cells used for computing integration vector weights, potentially making the integration process less sensitive to noise or variability in the data.

# Scale data, run PCA and UMAP and visualize integrated data
sobj_ys_filtered_integrated <- ScaleData(object = sobj_ys_filtered_integrated)
sobj_ys_filtered_integrated <- RunPCA(object = sobj_ys_filtered_integrated)
sobj_ys_filtered_integrated <- RunUMAP(object = sobj_ys_filtered_integrated, dims = 1:50)


p3 <- DimPlot(sobj_ys_filtered_integrated, reduction = 'umap', group.by = 'fetal.ids')
A1|p3

p4 <-DimPlot(sobj_ys_filtered_integrated, reduction = 'umap', group.by = 'LVL3')
A3 <- DimPlot(sobj_ys_filtered, reduction = 'umap', group.by = 'LVL3', label = T)

p4|A3

#making a 5had file for RNA Velocity

SaveH5Seurat(sobj_ys_filtered_integrated, filename = "sobj_ys_integrated.h5Seurat")
Convert("sobj_ys_integrated.h5Seurat", dest = "h5ad") 


#####MARIANA ADD ON's FOR PSEUDOTIME TRAJECTORY ANALYSIS######

#need to set assay in order to not run into issues later with monocle3
sobj_ys_filtered_integrated@active.assay='RNA'

# Join the layers of your Seurat object
sobj_ys_filtered_integrated <- JoinLayers(object = sobj_ys_filtered_integrated)

#converting seurat object to cell_data_set object for monocle3
cds <- SeuratWrappers::as.cell_data_set(sobj_ys_filtered_integrated)
cds

#get cell metadata
colData(cds)

#feature/gene metadata
fData(cds)
#has no columns 
rownames(fData(cds)) [1:10] #these are the genes

fData(cds)$gene_short_name <- rownames(fData(cds)) 

#to get counts

counts(cds)

##CLUSTERING CELLS (from Seurats UMAP)
#assigning partitions - using named vector

recreate.partition <- c(rep(1,length(cds@colData@rownames))) 
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition) #names are all the cells and the values are all 1

cds@clusters$UMAP$partitions <- recreate.partition #all cells being assigned to 1 partition can now see it stored 

#assigning cluster information

list_cluster <- sobj_ys_filtered@active.ident
cds@clusters$UMAP$clusters <- list_cluster

#Assigning UMAP Coordinate - cell embeddings

# Check the dimensions of existing UMAP dimension
umap_dim <- length(cds@int_colData@listData$reducedDims$UMAP)

cds@int_colData@listData$reducedDims$UMAP <- umap_data

cds@int_colData@listData$reducedDims$UMAP <- sobj_ys_filtered_integrated@reductions$umap@cell.embeddings #to store values inside of UMAP Plot

# plot

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,#want to place only one label on the cells that are carrying that label
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "LVL3",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names


#Trajectroy Graphs

cds <- learn_graph(cds, use_partition = FALSE) #know that cells follow one trajectroy so we set as false

plot_cells(cds,
           color_cells_by = 'LVL3',
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

# Ordering Cells

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 5]))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(stage, monocle3_pseudotime, median), fill = LVL2)) +
  geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(LVL3, monocle3_pseudotime, median), fill = LVL2)) +
  geom_boxplot()

# Finding Genes that Change as a Fxn of Pseudotime 
deg_mast <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

#deg_mast %>% 
#  arrange(q_value) %>% 
#  filter(status == 'OK') %>% 
#  head()

#FeaturePlot(sobj_ys_filtered_integrated, features = c('IL3','RUNX1','FLT3','CD34','ITGB7'))
#FeaturePlot(sobj_ys_filtered_integrated, features = c('KIT','GATA1','MITF','TNF','BATF','GATA2'))
FeaturePlot(sobj_ys_filtered_integrated, features = c('RUNX1','CD34','GATA2','GATA1','MITF','BATF','KIT','IL3'))

# visualizing pseudotime in seurat
sobj_ys_filtered_integrated$pseudotime <- pseudotime(cds)
Idents(sobj_ys_filtered_integrated) <- sobj_ys_filtered_integrated$LVL3
FeaturePlot(sobj_ys_filtered_integrated, features = "pseudotime", label = T)


# plotting subsets

cds_subset <- choose_cells(cds)

cds_subset

cds_subset_pt_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=8)
cds_subset_pt_res <- na.omit(cds_subset_pt_res)
cds_subset_pt_res <- cds_subset_pt_res[cds_subset_pt_res$p_value < 0.05 & cds_subset_pt_res$status == "OK", ]
cds_subset_pt_res

cds_subset_pt_res[order(-cds_subset_pt_res$morans_test_statistic),] 

plot_cells(cds_subset, genes=c('KIT','ITGB7','CD34'),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

cds_subset_subset <- cds_subset[rowData(cds_subset)$gene_short_name %in% c('GATA1','RUNX1','KIT','CD34','ITGB7')]

plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by="stage",
                         min_expr=0.5)

diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Media")
sig_genes <- subset(diff_test_res, qval < 0.1)

