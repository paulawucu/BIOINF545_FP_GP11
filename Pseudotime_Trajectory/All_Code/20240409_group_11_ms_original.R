##### LOTS OF INSTALLATION #####

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
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github('cole-trapnell-lab/monocle3')

library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(monocle3)
library(SeuratWrappers)

##### YS READ DATA #####

markers <- read.delim('ABC_Marker.txt', header = T) # gene metadata
metadata <- read.delim('metadata_ys.csv', header = T) # cell metadata
expr <- read.delim('matrix.mtx', header = T, sep = ',') # expression matrix

##### GENERATING SEURAT OBJECT #####
expr.t <- t(expr)
seu.obj <- CreateSeuratObject(counts = expr.t)
View(seu.obj@meta.data)
seu.obj@meta.data <- merge(seu.obj@meta.data, metadata, by.x = 'row.names', by.y = 'cell_id')
View(seu.obj@meta.data)
seu.obj@meta.data <- seu.obj@meta.data %>% 
  column_to_rownames(var = 'Row.names')
seu.obj$mitopercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')
seu.obj.filtered <- subset(seu.obj, subset = nCount_RNA > 800 &
                             nFeature_RNA > 500 &
                             mitopercent < 10)


####SUBSET SEURAT####

unique(seu.obj.filtered@meta.data$population)

Idents(seu.obj.filtered) <- seu.obj.filtered$population
b.seu <- subset(seu.obj.filtered, idents = "b")
b.seu
unique(b.seu@meta.data$redefined_cluster)

##### pre-processing using seurat #####
b.seu <- NormalizeData(b.seu)
b.seu <- FindVariableFeatures(b.seu)
b.seu <- ScaleData(b.seu)
b.seu <- RunPCA(b.seu)
b.seu <- FindNeighbors(b.seu, dims = 1:30)
b.seu <- FindClusters(b.seu, resolution = 0.9)
b.seu <- RunUMAP(b.seu, dims = 1:30, n.neighbors = 50)


#####MARIANA ADD ON's FOR PSEUDOTIME TRAJECTORY ANALYSIS######

#need to set assay in order to not run into issues later with monocle3
sobj_ys@active.assay='RNA'

#converting seurat object to cell_data_set object for monocle3
cds <- SeuratWrappers::as.cell_data_set(sobj_ys)
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

list_cluster <- sobj_ys@active.ident
cds@clusters$UMAP$clusters <- list_cluster

#Assigning UMAP Coordinate - cell embeddings

cds


