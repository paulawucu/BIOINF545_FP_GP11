#####ADDING LIBRARIES####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
install.packages("Seurat")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
install.packages("Matrix")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("devtools")
install.packages('dplyr')
install.packages("R.utils")
install.packages("harmony", dependencies = TRUE)
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages("remotes")
remotes::install_github("satijalab/Azimuth")

library(Seurat)
library(harmony)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(gridExtra)
library(Azimuth)
library(patchwork)

#####UPLOADING DATA#####
raw_data_ys <- Read10X(data.dir = "matrix_files_ys")

metadata_ys <- read.csv("metadata_ys.csv")

rownames(metadata_ys) <- metadata_ys$X

map_and_filter <- function(path_feature, raw_data_frame, metadata_frame){
  features <- read.delim(path_feature, header = FALSE, stringsAsFactors = FALSE)
  gene_ids <- features$V1
  gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  
  features$V2 <- gene_symbols
  features$V1 <- gene_symbols
  
  mapped_idx <- !is.na(gene_symbols)
  filtered_raw <- raw_data_frame[mapped_idx, ]
  filtered_features <- features[mapped_idx, ]
  
  rownames(filtered_raw) <- filtered_features$V1
  
  if(anyDuplicated(rownames(filtered_raw)) > 0){
    dup_ids <- which(duplicated(rownames(filtered_raw)))
    for(i in seq_along(dup_ids)){ old_name <- rownames(filtered_raw)[dup_ids[i]]
    new_name <- paste0(old_name, "_dup", i)
    rownames(filtered_raw)[dup_ids[i]] <- new_name
    }
  }
  
#####MAKING SEURAT OBJECT#####
  seuratObject <- CreateSeuratObject(counts = filtered_raw, meta.data = metadata_frame)
  return(seuratObject)
}

sobj_ys <- CreateSeuratObject(counts = raw_data_ys, meta.data = metadata_ys)

view(sobj_ys@meta.data)

sobj_ys$mitopercent <- PercentageFeatureSet(sobj_ys, pattern = "^MT-")

sobj_ys_filtered <- subset(sobj_ys, subset = nCount_RNA > 800 &
                             nFeature_RNA > 500 &
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

A1 <- DimPlot(sobj_ys_filtered, reduction = 'umap', group.by = 'LVL2', label = T)
  
A2 <- DimPlot(sobj_ys_filtered, reduction = 'umap', group.by = 'seurat_clusters', label = T)

A1|A2

#####Integration#####

##Integration Step using Harmony
#sobj_ys_integrated <- IntegrateLayers(object = sobj_ys_filtered, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",)
#obj_ys_integrated <- FindNeighbors(sobj_ys_integrated, reduction = "harmony", dims = 1:8)
#sobj_ys_integrated <- FindClusters(sobj_ys_integrated, resolution = 0.3, cluster.name = "harmony_cluster")
#sobj_ys_integrated <- RunUMAP(sobj_ys_integrated, reduction = "harmony", dims = 1:8, reduction.name = "harmony.umap")



#####MARIANA ADD ON's FOR PSEUDOTIME TRAJECTORY ANALYSIS######

#need to set assay in order to not run into issues later with monocle3
sobj_ys_filtered@active.assay='RNA'

#converting seurat object to cell_data_set object for monocle3
cds <- SeuratWrappers::as.cell_data_set(sobj_ys_filtered)
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

cds@int_colData@listData$reducedDims$UMAP <- sobj_ys_filtered@reductions$umap@cell.embeddings #to store values inside of UMAP Plot

# plot

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,#want to place only one label on the cells that are carrying that label
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "LVL2",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan', 'pink')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names


#Trajectroy Graphs

cds <- learn_graph(cds, use_partition = FALSE) #know that cells follow one trajectroy so we set as false

plot_cells(cds,
           color_cells_by = 'component',
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

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(component, monocle3_pseudotime, median), fill = LVL2)) +
  geom_boxplot()


# Finding Genes that Change as a Fxn of Pseudotime 
deg_mast <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_mast %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

#FeaturePlot(sobj_ys_filtered, features = c('E2F2', 'STMN1', 'CD52'))
FeaturePlot(sobj_ys_filtered, features = c('PTGR3', 'AC084033.3', 'H3-3A','H3-3B','COX1','ATP6'))
FeaturePlot(sobj_ys_filtered, features = c('GATA1','RUNX1','FLT3','CD34','ITGB7'))
FeaturePlot(sobj_ys_filtered, features = c('KIT','GATA1','MITF','TNF','BATF','IL33'))

# visualizing pseudotime in seurat

sobj_ys_filtered$pseudotime <- pseudotime(cds)
Idents(sobj_ys_filtered) <- sobj_ys_filtered$component
FeaturePlot(sobj_ys_filtered, features = "pseudotime", label = T)
