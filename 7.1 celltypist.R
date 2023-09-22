#我尝试用R来执行python
library('Matrix')
library('reticulate')
library('dplyr')
use_condaenv(condaenv = "py3.11", required = TRUE)#我自己的环境，你用的时候自己先在conda装一个，后面几个包提前在虚拟环境下用conda装好
mpl = import("matplotlib")
mpl$use("Agg")
scanpy = import("scanpy")
celltypist = import("celltypist")
pandas <- import("pandas")
numpy = import("numpy")
celltypist$models$download_models(force_update = F) 
counts_matrix <- Seurat_object[["RNA"]]@counts
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(counts_matrix))),
                       obs = pandas$DataFrame(Seurat_object@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(counts_matrix),
                                                         row.names = rownames(counts_matrix)))
)
#修改为可写，因为这个scanpy包好像只能只读
reticulate::py_run_string("
import numpy as np
adata.X = np.array(adata.X)
")
py_run_string("adata.X.setflags(write=True)")
# 将R中的adata对象分配给Python
py$adata <- adata
# 在Python中复制数据并进行归一化
reticulate::py_run_string("
import scanpy as sc
# 复制数据
adata_copy = sc.AnnData(X=adata.X.copy(), obs=adata.obs, var=adata.var)
# 在复制的数据上进行归一化
sc.pp.normalize_total(adata_copy, target_sum=1e4)
sc.pp.log1p(adata_copy)
")


scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
predictions = celltypist$annotate(adata, model = 'Human_AdultAged_Hippocampus.pkl', majority_voting = T)
predictions$predicted_labels %>% head()
Seurat_object = AddMetaData(Seurat_object, predictions$predicted_labels)
#在我不懈努力的debug下终于能运行了
