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
library(limma)

#####UPLOADING DATA#####
raw_data_ys <- Read10X(data.dir = "matrix_files_ys")

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

#sobj_ys_filtered@active.assay='RNA'

#class(sobj_ys_filtered)
#str(sobj_ys_filtered)


#####Removing Batch Effect#####

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

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)


#finding anchors 
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features,
                                  reduction = "rpca",
                                  dims = 1:3)
# integrate data
sobj_ys_integrated <- IntegrateData(anchorset = anchors, dims=1:3)
sobj_ys_integrated <- ScaleData(integrated)
sobj_ys_integrated <- RunPCA(integrated)
sobj_ys_integrated <- RunUMAP(integrated, dims = 1:10, reduction.name = "UMAP")
sobj_ys_integrated <- FindNeighbors(integrated, dims = 1:10)
spbj_ys_integrated <- FindClusters(integrated)
#DimPlot(integrated, group.by = c("orig.ident", "ident"))



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
                            color_cells_by = "fetal.ids",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan', 'pink','purple')) +
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

FeaturePlot(sobj_ys_filtered, features = c('GATA1','RUNX1','FLT3','CD34','ITGB7'))
FeaturePlot(sobj_ys_filtered, features = c('KIT','GATA1','MITF','TNF','BATF','IL33'))

# visualizing pseudotime in seurat

sobj_ys_filtered$pseudotime <- pseudotime(cds)
Idents(sobj_ys_filtered) <- sobj_ys_filtered$stage
FeaturePlot(sobj_ys_filtered, features = "pseudotime", label = T)


# plotting subsets

cds_subset <- choose_cells(cds)

cds_subset

cds_subset_pt_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=8)
cds_subset_pt_res <- na.omit(cds_subset_pt_res)
cds_subset_pt_res <- cds_subset_pt_res[cds_subset_pt_res$p_value < 0.05 & cds_subset_pt_res$status == "OK", ]
cds_subset_pt_res

cds_subset_pt_res[order(-cds_subset_pt_res$morans_test_statistic),] 

plot_cells(cds_subset, genes=c('GATA1','RUNX1','FLT3','CD34','ITGB7'),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

cds_subset_subset <- cds_subset[rowData(cds_subset)$gene_short_name %in% c('GATA1','RUNX1','FLT3','CD34','ITGB7')]

plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by="stage",
                         min_expr=0.5)


