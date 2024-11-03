library(GeneNMF)
library(Seurat)
library(tidyverse)
library(fgsea)
library(msigdbr)
setwd("/home/lfZhang/sc")
res_home <- "/home/lfZhang/sc/"

In <- readRDS("In_NMF.rds")
#geneNMF@reductions
geneNMF <- subset(In, subset = group == "MDD")
#geneNMF <- runNMF(geneNMF, assay = "SCT", k=30)
geneNMF <- SplitObject(geneNMF, split.by = "indiv")
#########################################
geneNMF.programs <- multiNMF(geneNMF, 
                             assay="SCT", slot="data", 
                             k=2:6, L1=c(0,0), 
                             do_centering=TRUE, 
                             nfeatures = 2000)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nprograms=10,
                                        max.genes=50,
                                        hclust.method="ward.D2",
                                        min.confidence=0.3)
geneNMF.metaprograms$metaprograms.metrics
ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,0.8))
ggsave(filename = "In_metaprograms.png", plot = ph,path = paste0(res_home, "Figure/"),  dpi = 400,width = 12, height = 12, limitsize = F)
t(as.data.frame(lapply(geneNMF.metaprograms$metaprograms.genes, head)))
#k=2:6,nfeatures = 2000,max.genes=50,
#     sampleCoverage  silhouette meanJaccard numberGenes numberPrograms
#MP1       1.0000000  0.24953155       0.396          50            118
#MP2       1.0000000  0.17261501       0.417          50             88
#MP3       0.9666667  0.18411659       0.289          50            101
#MP4       0.9333333  0.60175215       0.698          50             29
#MP5       0.9333333  0.12297298       0.329          50             84
#MP6       0.8333333  0.28015470       0.456          50             45
#MP7       0.6666667  0.04150070       0.150          42             40
#MP8       0.6666667 -0.01225382       0.065          10             44
#MP9       0.6333333  0.11400913       0.323          50             32
#MP10      0.3666667  0.30988374       0.386          50             19

universe <- rownames(In)
#saveRDS(universe, paste0(res_home,"universe_In.rds"))
#saveRDS(geneNMF.metaprograms,("geneNMF.metaprograms_In.rds"))
#saveRDS(geneNMF.programs,("geneNMF.programs_In.rds"))

universe <- readRDS("universe_In.rds")
geneNMF.metaprograms <- readRDS("geneNMF.metaprograms_In.rds")
geneNMF.programs <- readRDS("geneNMF.programs_In.rds")
##########C2############

program <- geneNMF.metaprograms$metaprograms.genes$MP2
library(dplyr)

top_p <- runGSEA(program, universe = universe, category = "C5") %>% 
  filter(padj < 0.05) %>%
  arrange(desc(overlap))
top_p <- top_p %>%
  mutate(pathway = sub("^[^_]*_", "", pathway)) %>%
  mutate(pathway = gsub("_", " ", pathway))
top_p$pathway <- str_wrap(top_p$pathway, width = 20) 
top_p <- top_p %>%
  top_n(n = 20, wt = overlap) 
########################

#########C5###########
top_p$group <- ifelse(grepl("GOBP", top_p$pathway), "BP",
                     ifelse(grepl("GOCC", top_p$pathway), "CC",
                            ifelse(grepl("GOMF", top_p$pathway), "MF", "HP")))

library(stringr)

top_p <- top_p %>%
  filter(group != "HP") %>%
  group_by(group) %>%
  top_n(n = 5, wt = overlap) %>%
  ungroup() %>% # 解除分组，以便后续操作不会受到分组影响
  mutate(pathway = sub("^[^_]*_", "", pathway))
top_p$pathway <- gsub("_", " ", top_p$pathway)
top_p$pathway <- str_wrap(top_p$pathway, width = 20) 
######################

##########figure############
# 加载ggplot2包
library(ggplot2)


ggplot(top_p, aes(x = reorder(pathway, overlap), y = overlap / geneNMF.metaprograms$metaprograms.metrics$numberGenes[4], fill = -log10(pval))) +
  geom_col(alpha = 0.7) +
  scale_fill_gradientn(colors = rev(viridis(10, option="mako"))) +
  theme_minimal() +
  coord_flip() +
  labs(
    fill = "-log10(p-value)"
  ) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(size = 15, face = "bold", angle = 90),
    strip.text = element_text(size = 15, face = "bold")
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1, strip.position = "right")
ggsave(paste0(res_home, "Figure/In_GO_MP6.png"),  dpi = 400,width = 10, height = 12, limitsize = F)


library(UCell)
mp.genes <- geneNMF.metaprograms$metaprograms.genes
In <- AddModuleScore_UCell(In, features = mp.genes, assay="SCT", ncores=4, name = "")
VlnPlot(In, features=names(mp.genes), group.by = "Types",
        pt.size = 0)
matrix <- In@meta.data[,names(mp.genes)]
dimred <- as.matrix(matrix)

colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
In@reductions[["MPsignatures"]] <- new("DimReduc",
                                         cell.embeddings = dimred,
                                         assay.used = "SCT",
                                         key = "MP_",
                                         global = FALSE)

In <- RunUMAP(In, reduction="MPsignatures", dims=1:length(In@reductions[["MPsignatures"]]),
               metric = "euclidean", reduction.name = "umap_MP")
library(viridis)
FeaturePlot(In, features = names(mp.genes), reduction = "umap", ncol=5) &
  scale_color_gradientn(colors = viridis(100, option="A")) &  
  theme(
    aspect.ratio = 1, 
    axis.title = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    axis.line = element_blank()
  )
ggsave(paste0(res_home, "Figure/In_MP.png"),  dpi = 400,width = 12, height = 4, limitsize = F)
#saveRDS(In, paste0(res_home,"In_NMF.rds"))
DimPlot(In, group.by = "Types", label = T, repel = T)




Ex <- readRDS("Ex_NMF.rds")
#geneNMF@reductions
geneNMF <- subset(Ex, subset = group == "MDD")
#geneNMF <- runNMF(geneNMF, assay = "SCT", k=30)
geneNMF <- SplitObject(geneNMF, split.by = "indiv")
#########################################
geneNMF.programs <- multiNMF(geneNMF, 
                             assay="SCT", slot="data", 
                             k=2:5, L1=c(0,0), 
                             do_centering=TRUE, 
                             nfeatures = 2000)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nprograms=10,
                                        max.genes=50,
                                        hclust.method="ward.D2",
                                        min.confidence=0.3)
geneNMF.metaprograms$metaprograms.metrics
ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,0.8))
#     sampleCoverage  silhouette meanJaccard numberGenes numberPrograms
#MP1       0.9333333  0.12257029       0.352          50             55
#MP2       0.9000000  0.17318330       0.434          50             32
#MP3       0.8666667  0.61968336       0.768          50             79
#MP4       0.8666667  0.11125106       0.290          50             51
#MP5       0.8333333  0.21987826       0.450          47             54
#MP6       0.7666667  0.80796928       0.893          50             28
#MP7       0.7666667 -0.07918276       0.077          14             51
#MP8       0.7000000  0.07406500       0.499          50             34
#MP9       0.4000000  0.04714029       0.421          50             21
#MP10      0.3000000  0.17897638       0.329          50             15

ggsave(filename = "Ex_metaprograms.png", plot = ph,path = paste0(res_home, "Figure/"),  dpi = 400,width = 12, height = 12, limitsize = F)
t(as.data.frame(lapply(geneNMF.metaprograms$metaprograms.genes, head)))
universe <- rownames(Ex)
#saveRDS(universe, paste0(res_home,"universe_Ex.rds"))
#saveRDS(geneNMF.metaprograms,("geneNMF.metaprograms_Ex.rds"))
#saveRDS(geneNMF.programs,("geneNMF.programs_Ex.rds"))

universe <- readRDS("universe_Ex.rds")
geneNMF.metaprograms <- readRDS("geneNMF.metaprograms_Ex.rds")
geneNMF.programs <- readRDS("geneNMF.programs_Ex.rds")
program <- geneNMF.metaprograms$metaprograms.genes$MP3
library(dplyr)

top_p <- runGSEA(program, universe = universe, category = "C2") %>% 
  filter(padj < 0.05) %>%
  arrange(desc(overlap)) #%>%
#  head(10)
#############C2###########

#########C5###########
top_p$group <- ifelse(grepl("GOBP", top_p$pathway), "BP",
                     ifelse(grepl("GOCC", top_p$pathway), "CC",
                            ifelse(grepl("GOMF", top_p$pathway), "MF", "HP")))


  

######################

##########figure############
# 加载ggplot2包
library(ggplot2)



########################

#########C5###########
top_p$group <- ifelse(grepl("GOBP", top_p$pathway), "BP",
                     ifelse(grepl("GOCC", top_p$pathway), "CC",
                            ifelse(grepl("GOMF", top_p$pathway), "MF", "HP")))

library(stringr)

top_p <- top_p %>%
  filter(group != "HP") %>%
  group_by(group) %>%
  top_n(n = 5, wt = overlap) %>%
  ungroup() %>% # 解除分组，以便后续操作不会受到分组影响
  mutate(pathway = sub("^[^_]*_", "", pathway))
top_p$pathway <- gsub("_", " ", top_p$pathway)
top_p$pathway <- str_wrap(top_p$pathway, width = 20) 
######################

##########figure############
# 加载ggplot2包
library(ggplot2)

ggplot(top_p, aes(x = reorder(pathway, overlap), y = overlap / geneNMF.metaprograms$metaprograms.metrics$numberGenes[4], fill = -log10(pval))) +
  geom_col(alpha = 0.7) +
  scale_fill_gradientn(colors = rev(viridis(10, option="rocket"))) +
  theme_minimal() +
  coord_flip() +
  labs(
    fill = "-log10(p-value)"
  ) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(size = 15, face = "bold", angle = 90),
    strip.text = element_text(size = 15, face = "bold")
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1, strip.position = "right")
ggsave(paste0(res_home, "Figure/Ex_GO_MP3.png"),  dpi = 400,width = 10, height = 12, limitsize = F)

library(UCell)
mp.genes <- geneNMF.metaprograms$metaprograms.genes
Ex <- AddModuleScore_UCell(Ex, features = mp.genes, assay="SCT", ncores=4, name = "")
VlnPlot(Ex, features=names(mp.genes), group.by = "subtypes",
        pt.size = 0)

library(viridis)
FeaturePlot(Ex, features = names(mp.genes), reduction = "umap", ncol=5) &
  scale_color_gradientn(colors = viridis(100, option="A")) &
    theme(
    aspect.ratio = 1, 
    axis.title = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    axis.line = element_blank()
  )
ggsave(paste0(res_home, "Figure/Ex_MP.png"),  dpi = 400,width = 12, height = 4, limitsize = F)
#saveRDS(Ex, paste0(res_home,"Ex_NMF.rds"))
#saveRDS(In, paste0(res_home,"In_NMF.rds"))


##################################
program <- geneNMF.metaprograms$metaprograms.genes$MP6
top_p <- runGSEA(program, universe = universe, category = "C2") %>% 
  filter(padj < 0.05) %>%
  arrange(desc(overlap))
top_p <- top_p %>%
  mutate(pathway = sub("^[^_]*_", "", pathway)) %>%
  mutate(pathway = gsub("_", " ", pathway))
top_p$pathway <- str_wrap(top_p$pathway, width = 20) 
top_p <- top_p %>%
  top_n(n = 20, wt = overlap)

library(stringr)
ggplot(top_p, aes(x = reorder(pathway, overlap), y = overlap / geneNMF.metaprograms$metaprograms.metrics$numberGenes[4], fill = -log10(pval))) +
  geom_col(alpha = 0.7) +
  scale_fill_gradientn(colors = rev(viridis(10, option="rocket"))) +
  theme_minimal() +
  coord_flip() +
  labs(
    fill = "-log10(p-value)"
  ) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(size = 15, face = "bold", angle = 90),
    strip.text = element_text(size = 15, face = "bold")
  )

ggsave(paste0(res_home, "Figure/Ex_KEGG_MP6.png"),  dpi = 400,width = 8, height = 18, limitsize = F)

program <- geneNMF.metaprograms$metaprograms.genes$MP6
top_p <- runGSEA(program, universe = universe, category = "C2") %>% 
  filter(padj < 0.05) %>%
  arrange(desc(overlap))
top_p <- top_p %>%
  mutate(pathway = sub("^[^_]*_", "", pathway)) %>%
  mutate(pathway = gsub("_", " ", pathway))
top_p$pathway <- str_wrap(top_p$pathway, width = 20) 
top_p <- top_p %>%
  top_n(n = 20, wt = overlap)

library(stringr)
ggplot(top_p, aes(x = reorder(pathway, overlap), y = overlap / geneNMF.metaprograms$metaprograms.metrics$numberGenes[4], fill = -log10(pval))) +
  geom_col(alpha = 0.7) +
  scale_fill_gradientn(colors = rev(viridis(10, option="mako"))) +
  theme_minimal() +
  coord_flip() +
  labs(
    fill = "-log10(p-value)"
  ) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(size = 15, face = "bold", angle = 90),
    strip.text = element_text(size = 15, face = "bold")
  )

ggsave(paste0(res_home, "Figure/In_KEGG_MP6.png"),  dpi = 400,width = 8, height = 18, limitsize = F)
