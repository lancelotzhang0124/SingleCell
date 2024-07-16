#To excute Python scripts in your R environment, please use package 'reticulate'
library('Matrix')
library('reticulate')
library('dplyr')
use_condaenv(condaenv = "py3.11", required = TRUE)
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
reticulate::py_run_string("
import numpy as np
adata.X = np.array(adata.X)
")
py_run_string("adata.X.setflags(write=True)")
py$adata <- adata
reticulate::py_run_string("
import scanpy as sc
adata_copy = sc.AnnData(X=adata.X.copy(), obs=adata.obs, var=adata.var)
sc.pp.normalize_total(adata_copy, target_sum=1e4)
sc.pp.log1p(adata_copy)
")


scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
predictions = celltypist$annotate(adata, model = 'Human_AdultAged_Hippocampus.pkl', majority_voting = T)
predictions$predicted_labels %>% head()
Seurat_object = AddMetaData(Seurat_object, predictions$predicted_labels)
