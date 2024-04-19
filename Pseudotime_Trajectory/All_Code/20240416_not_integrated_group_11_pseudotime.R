#####ADDING LIBRARIES####
# Install necessary packages if not already installed
packages_to_install <- c(
  "BiocManager", "Seurat", "Matrix", "ggplot2", "dplyr", "devtools", "R.utils", 
  "harmony", "remotes", "SeuratDisk","SeuratDisk","BioManager","monocle3"
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

BiocManager::install(version = "3.18")
BiocManager::install(c(
  'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'lme4', 
  'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 
  'HDF5Array', 'terra', 'ggrastr', 'AnnotationDbi', 'org.Hs.eg.db', 'DDRTree', 
  'pheatmap', 'DESeq2'
))

# Install packages from GitHub
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github("cole-trapnell-lab/DDRTree", ref="simple-ppt-like")
devtools::install_github("cole-trapnell-lab/L1-graph")
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
remotes::install_github("satijalab/Azimuth")


# Load libraries
# Data Analysis Libraries
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratWrappers)
library(monocle3)
library(DESeq2)
library(limma)
library(harmony)
library(ggplot2)
library(pheatmap)

# Bioinformatics Libraries
library(BiocManager)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Azimuth)

# Data Manipulation Libraries
library(Matrix)
library(tidyverse)
library(gridExtra)
library(dplyr)
library(magrittr)

#####UPLOADING DATA#####
getwd()
setwd("/Users/marianasierra/Desktop/Group_11/")

raw_data_ys <- Read10X(data.dir = "matrix_files_ys_meyloid_progen_cells")
metadata_ys <- read.csv("metadata_ys.csv")
rownames(metadata_ys) <- metadata_ys$X

#####MAKING SEURAT OBJECT#####

sobj_ys <- CreateSeuratObject(counts = raw_data_ys, meta.data = metadata_ys)
view(sobj_ys@meta.data)

#filter data for mitopercentage, nfeature RNA, and nCount
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
# Compute PCA verbose = FALSE controls whether the function displays progress and diagnostic messages during the principal component analysis (PCA) computation. 
sobj_ys_filtered <- RunPCA(sobj_ys_filtered, verbose = FALSE)
#Finding Neighboors
sobj_ys_filtered <- FindNeighbors(sobj_ys_filtered,dims = 1:30)
#clustering
sobj_ys_filtered <- FindClusters(sobj_ys_filtered, resolution = .9)
# Now run UMAP
sobj_ys_filtered <- RunUMAP(sobj_ys_filtered, dims = 1:30, n.neighbors = 50)

view(sobj_ys_filtered@meta.data)

#####Make Ratio Graphs and UMAPs to Check Data#####
# Get the sample and cell type information from metadata
sample_celltype <- data.frame(Sample = sobj_ys_filtered@meta.data$stage,
                              CellType = sobj_ys_filtered$LVL3)

# Calculate the frequencies of each cell type within each sample
celltype_table <- table(sample_celltype$Sample, sample_celltype$CellType)

# Calculate the ratios
celltype_ratios <- prop.table(celltype_table, margin = 1) * 100  # Convert to percentages

# Transpose the celltype_ratios table
celltype_ratios_transposed <- t(celltype_ratios)

# Define distinct colors for each cell type manually: When I used rainbow it didn't seem to have enough color to asign so I had to do it mannually
num_celltypes <- ncol(celltype_ratios_transposed)
colors <- c("red", "blue", "limegreen", "orange", "purple", "cyan", "magenta", "yellow", "lightgreen", "deeppink")

# Plot stacked bar graph
barplot(celltype_ratios_transposed, beside = FALSE,
        xlab = "Samples",
        ylab = "Percentage",
        main = "Percentage of Cell Types in Each Sample",
        col = colors,
        ylim = c(0, 100),  # Set y-axis limit from 0 to 100
        las = 2)  # Rotate x-axis labels vertically

# Add legend: I'm a weirdo and I prefer to have make the figure in powerpoint so I leave teh key off until then usually
#legend("topright", legend = rownames(celltype_ratios_transposed), fill = colors, title = "Cell Type",
#       xpd = TRUE,
#       text.width = 2)

#UMAPS
A1 <- DimPlot(sobj_ys_filtered, reduction = 'umap', group.by = 'fetal.ids', label = T)
A2 <- DimPlot(sobj_ys_filtered, reduction = 'umap', group.by = 'seurat_clusters', label = T)
A1|A2

ElbowPlot(sobj_ys_filtered)


#####MARIANA ADD ON's FOR PSEUDOTIME TRAJECTORY ANALYSIS######

####PREPROCESSING####
#need to set assay in order to not run into issues later with monocle3
sobj_ys_filtered@active.assay='RNA'
# Join the layers of your Seurat object
#sobj_ys_filtered <- JoinLayers(object = sobj_ys_filtered)
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

# Check the dimensions of existing UMAP dimension
#cds@int_colData@listData$reducedDims$UMAP <- sobj_ys_filtered_integrated@reductions$umap@cell.embeddings

umap_dim <- length(cds@int_colData@listData$reducedDims$UMAP)
cds@int_colData@listData$reducedDims$UMAP <- umap_data
cds@int_colData@listData$reducedDims$UMAP <- sobj_ys_filtered@reductions$umap@cell.embeddings #to store values inside of UMAP Plot

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

######Trajectroy Graphs#####
cds <- learn_graph(cds, use_partition = FALSE) #know that cells follow one trajectroy so we set as false
plot_cells(cds,
           color_cells_by = 'LVL3',
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

#Ordering Cells
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 3]))
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
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(stage, monocle3_pseudotime, median), fill = stage)) +
  geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(stage, monocle3_pseudotime, median), fill = LVL2)) +
  geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(stage, monocle3_pseudotime, median), fill = LVL3)) +
  geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(LVL3, monocle3_pseudotime, median), fill = LVL3)) +
  geom_boxplot()
# Finding Genes that Change as a Fxn of Pseudotime 
deg_mast <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)
FeaturePlot(sobj_ys_filtered, features = c('RUNX1','CD34','GATA2','GATA1','MITF','BATF','KIT','IL3',"FCER1A"))
# visualizing pseudotime in seurat
sobj_ys_filtered$pseudotime <- pseudotime(cds)
Idents(sobj_ys_filtered) <- sobj_ys_filtered$LVL3
FeaturePlot(sobj_ys_filtered, features = "pseudotime", label = T)
#####SubSetting####
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
cds_subset_subset <- cds_subset[rowData(cds_subset)$gene_short_name %in% c('GATA1','RUNX1','KIT','CD34','ITGB7','FCER1A')]
plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by="stage",
                         min_expr=0.5)
plot_genes_in_pseudotime(cds_subset_subset,
                         color_cells_by="LVL3",
                         min_expr=0.5)


#####DEGS MAP MONOCLE 3/DESeq2#####

#expression matrix for cds
expression_matrix <- exprs(cds)
expression_matrix <- round(expression_matrix)

dds <- DESeqDataSetFromMatrix(countData = expression_matrix,
                              colData = sobj_ys_filtered@meta.data,
                              design = ~ pseudotime)
dds <- DESeq(dds)

dds$stage<- factor(dds$stage)

# Perform differential expression analysis to identify variable genes
res_var <- results(dds, contrast = c("stage","CS10","CS23"))

# Select top variable genes
top_variable_genes <- rownames(head(res_var[order(res_var$padj), ], 10))

# Extract expression data for the top variable genes
variable_genes_expr <- expression_matrix[top_variable_genes, ]
