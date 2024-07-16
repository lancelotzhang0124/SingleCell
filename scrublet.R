library('Matrix')
library('reticulate')
library('Seurat')
counts_matrix <- Seurat_object[["SCT"]]@counts
py$counts_matrix <- t(counts_matrix)
reticulate::py_run_string("
import numpy as np
import pandas as pd
import scrublet as scr
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
")
scrublet <- py$predicted_doublets
Seurat_object@meta.data$scrublet <- scrublet
str(scrublet)
true_count <- sum(scrublet == "TRUE", na.rm = TRUE)
saveRDS(Seurat_object,"Seurat_object_qc.rds")