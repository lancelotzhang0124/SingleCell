#############################################
library(Seurat)
library(GeneNMF)
library(ggplot2)
library(viridis)
seu <- readRDS("seu_NMF.rds")
.mp_cols  <- c(
                     "#c28948",
                     "#4bbbb5",
                     "#1c5294", 
                     "#d85bb3",
                     "#5d7fe5",
                     "#b4492e",
                     "#3aa176",
                     "#5d7b4e",
                     "#733a83",
                     "#ac5f72")
seu_nmf <- subset(seu, subset = group == "MDD")
seu_nmf <- SplitObject(seu_nmf, split.by = "indiv")

geneNMF.programs <- multiNMF(seu_nmf, 
                             assay="SCT", slot="data", 
                             k=2:4, L1=c(0, 0),
                             nfeatures = 2000)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        max.genes=100,
                                        hclust.method="ward.D2",
                                        nMP=10,
                                        metric = c("cosine", "jaccard"),
                                        specificity.weight = 5,
                                        weight.explained = 0.7,
                                        min.confidence=0.3)

print(geneNMF.metaprograms$metaprograms.metrics)

# function annotation
functional_names <- c(
    "_Neuronal signaling", 
    "_Immune response", 
    "_Cell proliferation", 
    "_Metabolism", 
    "_Synaptic transmission",
    "_Cell adhesion",
    "_Developmental processes",
    "_Inflammatory response",
    "_Signal transduction",
    "_Cytoskeleton organization"
)
cluster_to_function <- setNames(functional_names, 1:10)

geneNMF.metaprograms[["programs.clusters"]][!is.na(geneNMF.metaprograms[["programs.clusters"]])] <- 
    cluster_to_function[geneNMF.metaprograms[["programs.clusters"]][!is.na(geneNMF.metaprograms[["programs.clusters"]])]]
str(geneNMF.metaprograms)

# define colors
annotation_colors <- list(
    Metaprogram = setNames(
        .mp_cols,
        paste0("MP", functional_names)
    )
)
print(annotation_colors)

ph <- plotMetaPrograms(geneNMF.metaprograms, palette = colorRampPalette(c("#e2e7e9", "#0f2b5a"))(100), annotation_colors = annotation_colors, jaccard.cutoff = c(0.1,1))
