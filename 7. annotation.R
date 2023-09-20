
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
