if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("Seurat")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
install.packages("Matrix")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("devtools")
install.packages('dplyr')

library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(dplyr)

raw_data_human <- Read10X(data.dir = "matrix_files_humans")
raw_data_ys <- Read10X(data.dir = "matrix_files_ys")

metadata_human <- read.csv("metadata_human.csv")
metadata_ys <- read.csv("metadata_ys.csv")

rownames(metadata_human) <- metadata_human$X
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
  
  seuratObject <- CreateSeuratObject(counts = filtered_raw, meta.data = metadata_frame)
  return(seuratObject)
}

sobj_adult_noNA = map_and_filter('matrix_files_humans/features.tsv.gz', raw_data_humans, metadata_human)

sobj_human <- CreateSeuratObject(counts = raw_data_humans, meta.data = metadata_human)
sobj_ys <- CreateSeuratObject(counts = raw_data_ys, meta.data = metadata_ys)

saveRDS(sobj_adult_noNA, "sobj_adult.rds")
saveRDS(sobj_ys, "sobj_ys.rds")

sobj_adult_noNA$orig.ident <- "adult"
sobj_ys$orig.ident <- "yolk_sac"

all <- merge(x=sobj_adult_noNA, y = sobj_ys, add.cell.ids = c("adult", "yolk_sac"), project = "combine_clustering")

ages <- ifelse(
  is.na(all@meta.data$development_stage),
  0,
  as.numeric(str_extract(all@meta.data$development_stage, "\\d+"))  # Extract numeric age
)
all <- AddMetaData(object = all, metadata = ages, col.name = "age")

all[[]]

all <- NormalizeData(all)
all <- FindVariableFeatures(all)
all <- ScaleData(all, vars.to.regress = c("sex","age"))

all <- RunPCA(all)
all <- FindNeighbors(all, dims = 1:10, reduction = "pca")

all <- FindClusters(all, resolution = 0.3, cluster.name = "unintegrated_clusters")

all <- RunUMAP(all, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(all, reduction = "umap.unintegrated",group.by = "seurat_clusters", split.by = "orig.ident")
DimPlot(all, reduction = "umap.unintegrated", group.by = c("orig.ident"))
DimPlot(all, reduction = "umap.unintegrated", group.by = c("seurat_clusters","tissue_in_publication"), split.by = "orig.ident")

all_integrated <- IntegrateLayers(object = all, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",)
all_integrated <- FindNeighbors(all_integrated, reduction = "harmony", dims = 1:10)
all_integrated <- FindClusters(all_integrated, resolution = 0.25, cluster.name = "harmony_cluster")
all_integrated <- RunUMAP(all_integrated, reduction = "harmony", dims = 1:10, reduction.name = "harmony.umap")

DimPlot(all_integrated, reduction = "harmony.umap", group.by = c("orig.ident", "harmony_cluster"))
DimPlot(all_integrated, reduction = "harmony.umap", group.by = c("harmony_cluster","tissue_in_publication"), split.by = "orig.ident")
Dimplot(all_integrated, reduction = "harmony.umap", group.by = c("orig.ident", "harmony_cluster"))

saveRDS(all_integrated, file = "integrated_data.rds")
saveRDS(all, file = "unitegrated_integrated.rds")


      