###########################标准化#####################
Seurat_object <- NormalizeData(Seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
###########################高变基因#####################
Seurat_object <- FindVariableFeatures(Seurat_object, 
                                      selection.method = "vst", 
                                      nfeatures = 2000)
##########################归一化#####################
Seurat_object <- ScaleData(Seurat_object, features = rownames(Seurat_object))
#以上每一步都需要大量运算能力，个人电脑慎运行
