library(msigdbr)
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
res_home <- "/home/lfzhang/SingleCell/"
scRNA_harmony = Ex

# First, let's look at the annotated UMAP plot
DimPlot(scRNA_harmony, reduction = "umap", group.by = "SCT_snn_res.0.3", label = TRUE, pt.size = 0.3, repel = TRUE)

# Check the number of cells in each cluster
scRNA_harmony$seurat_clusters <- scRNA_harmony$SCT_snn_res.0.3
table(scRNA_harmony$SCT_snn_res.0.3)

#### Differential Analysis ####
cell.markers <- Ex_0.3
Idents(scRNA_harmony) = "SCT_snn_res.0.3"
#1. Find markers for each cluster compared to all remaining cells
cell.markers <- FindAllMarkers(object = scRNA_harmony, 
                               only.pos = FALSE, # Report both upregulated and downregulated genes
                               test.use = "wilcox", # Default is wilcox test, see documentation for other options
                               slot = "data", # Note that the default is data, not counts
                               min.pct = 0.25, # Set the threshold for expression percentage, adjust higher if too many differential genes
                               logfc.threshold = 0.25 # Set the log2FC threshold, adjust higher if too many differential genes
)

# View the distribution of cells and clusters
table(scRNA_harmony@meta.data$MainType, scRNA_harmony@meta.data$seurat_clusters)

# 2. Differential analysis between two clusters
# pct.1 represents the proportion of single-cell samples expressing the gene in cluster 1, and pct.2 represents the proportion in cluster 2

Ex_degs = FindMarkers(scRNA_harmony, 
                      logfc.threshold = 0.25,
                      min.pct = 0.1, # Set a smaller threshold to find more differential genes
                      only.pos = FALSE,
                      ident.1 = "T cells", ident.2 = "AT2") %>% # ident 1 vs 2
  mutate(gene = rownames(.)) # Add row names as a new column for easier gene retrieval


#### GO/KEGG Enrichment ####
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(dplyr)
C0_degs_fil <- Ex_0.3 %>% filter(cluster == 0)

# Add corresponding ENTREZID for each gene
ids = bitr(C0_degs_fil$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db')

# Merge data, genes without ENTREZID in cluser3.markers will be filtered out
C0_degs_fil = merge(C0_degs_fil, ids, by.x = 'gene', by.y = 'SYMBOL')

# View data structure
head(C0_degs_fil)

# Sort genes by avg_log2FC in descending order
C0_degs_fil <- C0_degs_fil[order(C0_degs_fil$avg_log2FC, decreasing = TRUE), ]

# Generate gene list containing only ENTREZID names and avg_log2FC values
C0_degs_list <- as.numeric(C0_degs_fil$avg_log2FC)
names(C0_degs_list) <- C0_degs_fil$ENTREZID
head(C0_degs_list)

# enrichGO
# Select genes with larger differences
C0_de <- names(C0_degs_list)[abs(C0_degs_list) > 0.5]
head(C0_de)

# GO enrichment
C0_ego <- enrichGO(C0_de, OrgDb = "org.Hs.eg.db", ont = "BP", readable = TRUE)
head(C0_ego)

# Bubble plot
dotplot(C0_ego, showCategory = 10, title = "cluster0")
ggsave(filename = "C0_enrich_GO.png", path = paste0(res_home, "Figure/"), height = 5, width = 7)

# KEGG enrichment
C0_ekg <- enrichKEGG(gene = C0_de, organism = "hsa", pvalueCutoff = 0.05)
head(C0_ekg)

# Bubble plot
dotplot(C0_ekg, showCategory = 10, title = "cluster0 KEGG")
ggsave(filename = "C0_enrich_KEGG.png", path = paste0(res_home, "Figure/"), height = 4, width = 7)

C1_degs_fil <- Ex_0.3 %>% filter(cluster == 1)

# Add corresponding ENTREZID for each gene
ids = bitr(C1_degs_fil$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db')

# Merge data, genes without ENTREZID in cluser3.markers will be filtered out
C1_degs_fil = merge(C1_degs_fil, ids, by.x = 'gene', by.y = 'SYMBOL')

# View data structure
head(C1_degs_fil)

# Sort genes by avg_log2FC in descending order
C1_degs_fil <- C1_degs_fil[order(C1_degs_fil$avg_log2FC, decreasing = TRUE), ]

# Generate gene list containing only ENTREZID names and avg_log2FC values
C1_degs_list <- as.numeric(C1_degs_fil$avg_log2FC)
names(C1_degs_list) <- C1_degs_fil$ENTREZID
head(C1_degs_list)

results <- list()
enrich_results <- list()
i = 0

for (i in 0:12) {
  # Use tryCatch to handle possible errors
  tryCatch({
    # Filter data for specific cluster
    degs_fil <- Ex_0.3 %>% filter(cluster == i)
    
    # Check if degs_fil is empty
    if (nrow(degs_fil) == 0) stop("No data for cluster ", i)
    # If not cluster 0, perform replacement operation
    degs_fil$gene <- sub("\\..*", "", degs_fil$gene)
    
    # Add corresponding ENTREZID for each gene
    ids <- bitr(degs_fil$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db')
    
    # Merge data, genes without ENTREZID will be filtered out
    degs_fil <- merge(degs_fil, ids, by.x = 'gene', by.y = 'SYMBOL')
    
    # Select genes with larger differences
    degs_list <- setNames(as.numeric(degs_fil$avg_log2FC), degs_fil$ENTREZID)
    de <- names(degs_list)[abs(degs_list) > 0.5]
    
    # GO enrichment
    ego <- enrichGO(de, OrgDb = 'org.Hs.eg.db', ont = "BP", readable = TRUE)
    enrich_results[[paste0("C", i, "_ego")]] <- ego
    
    # KEGG enrichment
    ekg <- enrichKEGG(gene = de, organism = "hsa", pvalueCutoff = 0.05)
    enrich_results[[paste0("C", i, "_ekg")]] <- ekg
    
    # Plot and save bubble plot
    if (!is.null(ego)) {
      dotplot(ego, showCategory = 10, title = paste0("cluster", i, " GO"))
      ggsave(filename = paste0("C", i, "_enrich_GO.png"), path = paste0(res_home, "Figure/"), height = 5, width = 7)
    }
    
    if (!is.null(ekg)) {
      dotplot(ekg, showCategory = 10, title = paste0("cluster", i, " KEGG"))
      ggsave(filename = paste0("C", i, "_enrich_KEGG.png"), path = paste0(res_home, "Figure/"), height = 4, width = 7)
    }
  }, error = function(e) {
    cat("An error occurred in cluster", i, ": ", e$message, "\n")
  })
}

C7_degs_fil <- Ex_0.3 %>% filter(cluster == 7)

# Add corresponding ENTREZID for each gene
C7_degs_fil$gene <- sub("\\..*", "", rownames(C7_degs_fil))

ids = bitr(C7_degs_fil$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db')

# Merge data, genes without ENTREZID in cluser3.markers will be filtered out
C7_degs_fil = merge(C7_degs_fil, ids, by.x = 'gene', by.y = 'SYMBOL')

# View data structure
head(C7_degs_fil)

# Sort genes by avg_log2FC in descending order
C7_degs_fil <- C7_degs_fil[order(C7_degs_fil$avg_log2FC, decreasing = TRUE), ]

# Generate gene list containing only ENTREZID names and avg_log2FC values
C7_degs_list <- as.numeric(C7_degs_fil$avg_log2FC)
names(C7_degs_list) <- C7_degs_fil$ENTREZID
head(C7_degs_list
