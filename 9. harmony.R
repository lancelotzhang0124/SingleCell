packages <- c("Seurat","RColorBrewer","ggthemes","ggplot2","tidydr","ggalluvial","CellChat","harmony","viridis")
lapply(packages, function(i){library(i,character.only = T)})
res_home <- "/home/lfZhang/SingleCell/"

Seurat_object <- readRDS('Seurat_object.rds')
#######################data loading###########################
colnames(Seurat_object@meta.data)
dim(Seurat_object@meta.data)
table(Seurat_object@meta.data$seurat_clusters)
table(Seurat_object@meta.data$orig.ident)
table(Seurat_object@meta.data$predicted_labels)

.cluster_cols <- c("#ffc1dc", "#d76092", "#d03d33", "#ff3c28", "#ca0123", 
                   "#d11fa1", "#7f48c3", "#7b68ee", "#191970", "#1a3d85",
                   "#3d68c4", "#507dfe", "#0097ce", "#94cad0", "#4ab78d",
                   "#a1ff59", "#37864d", "#ffeebc", "#ffca15", "#f79d00", 
                   "#d07733", "#ff7832")

#先粗略的看一下做出来的图是否存在明显的批次效应
Idents(Seurat_object) <- Seurat_object$seurat_clusters
pdf(paste0(res_home, "Figure/Seurat_object_UMAP_seurat_clusters.pdf"),width = 8,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$seurat_clusters))),label = F,repel = T,label.size = 3,label.box = T) & theme_dr() & theme(panel.grid=element_blank())
dev.off()

Idents(Seurat_object) <- Seurat_object$predicted_labels
pdf(paste0(res_home, "Figure/Seurat_object_UMAP_predicted_labels.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$predicted_labels))),label = F,repel = T,label.size = 3,label.box = T) & theme_dr() & theme(panel.grid=element_blank())
dev.off()

###############batch correction###############
#需要鸿蒙
library(harmony)
library(devtools)

Seurat_object@meta.data$predicted_labels <- as.factor(Seurat_object@meta.data$predicted_labels)

Seurat_object_harmony <- RunHarmony(Seurat_object ,"predicted_labels", plot_convergence = TRUE)
harmony_embeddings<- Embeddings(Seurat_object_harmony, "harmony")
library(tidyverse)
Seurat_object_harmony <- Seurat_object_harmony %>% RunUMAP(reduction="harmony", dims = 1:20) %>% FindNeighbors(reduction = "harmony", dims = 1:20)
Seurat_object <- FindClusters(Seurat_object_harmony, resolution = seq(0.2, 0.5, by=0.1), algorithm = 1)
#这一步已经替换了原来的Seurat_object里的0.2~0.5分辨率的cluster
colnames(Seurat_object@meta.data)

saveRDS(Seurat_object_harmony, file = 'Seurat_object_harmony.rds')



Idents(Seurat_object) <- Seurat_object$seurat_clusters
pdf(paste0(res_home, "Figure/Seurat_object_UMAP_seurat_clusters.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$seurat_clusters))),label = F,repel = T,label.size = 3,label.box = T) & theme_dr() & theme(panel.grid=element_blank())
dev.off()

Idents(Seurat_object) <- Seurat_object$predicted_labels
pdf(paste0(res_home, "Figure/Seurat_object_UMAP_predicted_labels.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$seurat_clusters))),label = F,repel = T,label.size = 3,label.box = T) & theme_dr() & theme(panel.grid=element_blank())
dev.off()

Idents(Seurat_object) <- Seurat_object$RNA_snn_res.0.5
pdf(paste0(res_home, "Figure/Seurat_object_UMAP_RNA_snn_res.0.5.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$seurat_clusters))),label = F,repel = T,label.size = 3,label.box = T) & theme_dr() & theme(panel.grid=element_blank())
dev.off()

Idents(Seurat_object) <- Seurat_object$RNA_snn_res.0.4
pdf(paste0(res_home, "Figure/Seurat_object_UMAP_RNA_snn_res.0.4.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$seurat_clusters))),label = F,repel = T,label.size = 3,label.box = T) & theme_dr() & theme(panel.grid=element_blank())
dev.off()

Idents(Seurat_object) <- Seurat_object$RNA_snn_res.0.3
pdf(paste0(res_home, "Figure/Seurat_object_UMAP_RNA_snn_res.0.3.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$seurat_clusters))),label = F,repel = T,label.size = 3,label.box = T) & theme_dr() & theme(panel.grid=element_blank())
dev.off()

Idents(Seurat_object) <- Seurat_object$RNA_snn_res.0.2
pdf(paste0(res_home, "Figure/Seurat_object_UMAP_RNA_snn_res.0.2.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$seurat_clusters))),label = F, repel = T,label.size = 3,label.box = T) & theme_dr() & theme(panel.grid=element_blank())
dev.off()

Idents(Seurat_object) <- Seurat_object$RNA_snn_res.0.2
pdf(paste0(res_home, "Figure/Seurat_object_UMAP_RNA_snn_res.0.2_harmony.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'harmony',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$seurat_clusters))),label = T, repel = T,label.size = 3,label.box = T) & theme_dr() & theme(panel.grid=element_blank())
dev.off()
