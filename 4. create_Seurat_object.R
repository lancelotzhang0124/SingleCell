#####################构建 Seurat 对象##################### 
fixed_feature_names <- gsub("_", "-", feature.names)
rownames(mat) = fixed_feature_names
Seurat_object <- CreateSeuratObject(
  counts = mat, # 表达矩阵，可以为稀疏矩阵，也可以为普通矩阵
  min.cells = 3, # 去除在小于3个细胞中表达的基因
  min.features = 200) # 去除只有 200 个以下基因表达的细胞
gene_exp = Seurat_object@assays[["RNA"]]@counts
meta.data = Seurat_object@meta.data
Seurat_object[["percent.mt"]] <- PercentageFeatureSet(Seurat_object, pattern = "^MT")
Seurat_object <- subset(Seurat_object, subset = nCount_RNA > 1000 & 
                          nFeature_RNA < 5000 & #可以自行定义选择标准
                          percent.mt < 5 & 
                          nFeature_RNA > 600)
