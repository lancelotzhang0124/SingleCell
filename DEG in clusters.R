##############cluster间差异分析#####################

#cluster6.markers <- FindMarkers(Seurat_object, ident.1 = 6, min.pct = 0.25)
#找出区分cluster 5与cluster6的标记
#cluster5.markers <- FindMarkers(Seurat_object, ident.1 = 5, ident.2 = 6, min.pct = 0.25)
# 找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#FindAllMarkers函数用于找出每个聚类相对于其他所有细胞的标志性基因，并只报告正向的标志性基因（即在给定聚类中表达量高于其他细胞的基因）。
#接着，它使用dplyr的管道操作，对每个聚类筛选出平均log2折叠变化（avg_log2FC）最高的两个基因。
#VlnPlot(Seurat_object, features = c("USP11"))
top_genes_per_cluster <- all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top_genes_per_cluster

pdf("p6.pdf", width=10, height=20)
DimHeatmap(Seurat_object, dims = 1:15, cells = 500, balanced = TRUE) #?个PC 500个细胞
dev.off()
#未完成
