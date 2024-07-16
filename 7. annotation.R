
###################SingleR annotation celltypes###################
library(SingleR)
library(celldex)
#obtain annotation data base
hpca.se <- HumanPrimaryCellAtlasData()#create file .cache/R/ExperimentHub
hpca.se
Standard_matrix <- GetAssayData(Seurat_object, slot="data")
clusters=Seurat_object@meta.data$seurat_clusters
#annotation
pred.hpca.se <- SingleR(test = Standard_matrix, ref = hpca.se, labels = hpca.se$label.main,
                              method = "cluster", clusters = clusters,
                              assay.type.test = "logcounts",
                              assay.type.ref = "logcounts") 
cellType=data.frame(ClusterID=levels(Seurat_object@meta.data$seurat_clusters),
                    celltype=pred.hpca.se$labels)
Seurat_object@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),'celltype']
