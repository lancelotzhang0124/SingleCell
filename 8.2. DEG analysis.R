##################暴露间差异分析####################
Seurat_object <- SetIdent(Seurat_object, value = "group") #这里改变了ident了，切记！！！
#Seurat_object <- SetIdent(Seurat_object, value = "seurat_clusters")#还原ident
#寻找差异markers
diff_genes <- FindMarkers(Seurat_object, ident.1 = "Control", ident.2 = "Suicide")
top_genes <- diff_genes %>% 
  arrange(desc(abs(avg_log2FC))) %>%  
  head(10) 
gene_name <- rownames(top_genes)
print(gene_name)
