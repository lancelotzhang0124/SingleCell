
###################SingleR自动注释细胞###################
library(SingleR)
library(celldex)
#获取数据库
hpca.se <- HumanPrimaryCellAtlasData()#这里会自动创建一个文档.cache/R/ExperimentHub
hpca.se
#获取标准化矩阵
Standard_matrix <- GetAssayData(Seurat_object, slot="data")
#定义clusters变量
clusters=Seurat_object@meta.data$seurat_clusters
#注释细胞
pred.hpca.se <- SingleR(test = Standard_matrix, ref = hpca.se, labels = hpca.se$label.main,
                              method = "cluster", clusters = clusters,
                              assay.type.test = "logcounts",
                              assay.type.ref = "logcounts") 
cellType=data.frame(ClusterID=levels(Seurat_object@meta.data$seurat_clusters),
                    celltype=pred.hpca.se$labels)
Seurat_object@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),'celltype']


####################手动注释细胞####################
#由于我的数据已经标注好了细胞，并存在orig.ident中
# 1. 我们为每个聚类找到最常见的orig.ident
most_common_celltypes <- sapply(levels(Seurat_object$seurat_clusters), function(cluster) {
  cluster_cells <- subset(Seurat_object@meta.data, seurat_clusters == cluster)
  common_celltype <- names(sort(table(cluster_cells$orig.ident), decreasing = TRUE)[1])
  return(common_celltype)
})

# 2. 匹配聚类到细胞类型
Seurat_object$predicted_celltype <- Seurat_object$seurat_clusters
levels(Seurat_object$predicted_celltype) <- most_common_celltypes

# 输出聚类与细胞类型的对应关系
print(most_common_celltypes)

