################loading###############
library("dplyr")
library("cowplot")
library("stringr")
library("ggplot2")
library("ggrepel")
library("Seurat")
library('reticulate')
library("harmony")
library("RColorBrewer")
library("viridis")
setwd("/home/lfzhang/SingleCell/")
res_home <- "/home/lfzhang/SingleCell/"
.cluster_cols <- c("#92454e","#d8c9e3","#d4c380","#283b95","#385b7e")# Paul Signac 1

c("#dacc57","#c67d49","#985323","#14357b","#1b599c","#8bb6c2")# Vincet van Gogh
c("#f2d2cd","#e4ce84","#c3556f","#7d7db2","#65a99e")# Paul Signac 1
c("#d4c84b","#dd8f1c","#b36443","#72a751","#3c5b90")# Paul Signac 2
c("#d8c9e3","#d4c380","#92454e","#385b7e","#283b95")# Paul Signac 3
c("#92454e","#d8c9e3","#d4c380","#283b95","#385b7e")# Claude Monet 1
c("#b79f53","#9f9bbb","#ac5236","#5f2f33","#5d7b4e")# Claude Monet 2

##########multisession#####################
library("future")
library("future.apply")
plan(multisession)
options(future.globals.maxSize = 512 * 1024^3)
future::plan(multisession, workers=20)
future::plan(sequential)
#########expression matrix#####################
matrix_dir = "/home/lfzhang/SingleCell/GSEdata/"
barcode.path <- paste0(matrix_dir, "GSE213982_combined_counts_matrix_cells_columns.csv.gz")
features.path <- paste0(matrix_dir, "GSE213982_combined_counts_matrix_genes_rows.csv.gz")
matrix.path <- paste0(matrix_dir, "GSE213982_combined_counts_matrix.mtx")
mat <- readMM(file = matrix.path)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names <- barcode.names[-1, ]#deleted rownames
barcode.names <- sub("^\\d+,", "", barcode.names)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
feature.names <- feature.names[-1, ]
feature.names <- sub("^\\d+,", "", feature.names)
colnames(mat) = barcode.names
rownames(mat) = feature.names
mat
#####################create Seurat object#####################
fixed_feature_names <- gsub("_", "-", feature.names)
rownames(mat) = fixed_feature_names
Seurat_object <- CreateSeuratObject(counts = mat)
####################quality control######################
Seurat_object[["percent.hb"]] <- PercentageFeatureSet(Seurat_object, pattern = "^HB")
rb.genes <- rownames(Seurat_object)[grepl('^RP',rownames(Seurat_object))]
C<-GetAssayData(object = Seurat_object, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
Seurat_object <- AddMetaData(Seurat_object, percent.ribo, col.name = "percent.ribo")
p.rb <- VlnPlot(Seurat_object, features = c("percent.ribo"), pt.size = 0,log = TRUE, ncol = 2, raster = F)
ggsave(filename = "rb.png", path = paste0(res_home, "Figure/"), plot = p.rb,  dpi = 400, width = 15, height = 10, limitsize = F)
rm(C)
Seurat_object <- subset(Seurat_object, (subset = nCount_RNA > 350 & 
                          nCount_RNA < 25000 &
                          nFeature_RNA < 5000 &
                          percent.ribo < 30 &
                          nFeature_RNA > 200))
####################group info######################
#chemistry batch
group_values <- ifelse(str_detect(colnames(Seurat_object), "F18|F19|F20|F27|F28|F9|M1|M10|M11|M14|M17|M18|M23|M26|M28|M30|M32|M33|M34|M4|M5|M6|M8|F23|F24|F26|F29|F32|M12|M13|M15|M16|M19|M2|M20|M21|M22|M24|M24_2|M27|M29|M3|M31|M7|M9"
), "v2", 
                       ifelse(str_detect(colnames(Seurat_object), "F1|F2|F3|F4|F5|F6|F7|F8|F9|F11|F12|F14|F15|F16|F17|F18|F19|F20|F25|F27|F28|M1|M4|M5|M6|M8|M10|M11|M14|M17|M18|M23|M26|M28|M30|M32|M33|M34
"), "v3", NA))
group_values
Seurat_object$chem <- group_values
as.factor(Seurat_object$chem)
table(Seurat_object$chem)
#groups
group_values <- ifelse(str_detect(colnames(Seurat_object), "F7|F10|F13|F21|F22|F23|F24|F26|F29|F30|F31|F32|F33|F34|F35|F36|F37|F38|M2|M3|M7|M9|M12|M13|M15|M16|M19|M20|M21|M22|M24|M24_2|M27|M29|M31
"), "HC", 
                       ifelse(str_detect(colnames(Seurat_object), "F1|F2|F3|F4|F5|F6|F8|F9|F11|F12|F14|F15|F16|F17|F18|F19|F20|F25|F27|F28|M1|M4|M5|M6|M8|M10|M11|M14|M17|M18|M23|M26|M28|M30|M32|M33|M34
"), "MDD", NA))
group_values
Seurat_object$group <- group_values
Idents(Seurat_object) <- Seurat_object$group
p.group <- DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$group))), group.by = "group", label = F, repel = T, label.size = 3, label.box = TRUE, raster = F)
p.group
#individual batch
string <- colnames(Seurat_object)
indiv_name <- sub("\\..*$", "", string)
Seurat_object$indiv <- indiv_name
p.indiv <- DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$batch))), group.by = "batch", label = F, repel = T, label.size = 3, label.box = TRUE, raster = F)
p.indiv
#sex groups
sex_values <- ifelse(str_detect(colnames(Seurat_object), "F"), "Female", 
                     ifelse(str_detect(colnames(Seurat_object), "M"), "Male", NA))
sex_values
Seurat_object$sex <- sex_values
p.sex <- DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(length(table(Seurat_object$sex))), group.by = "sex", label = F, repel = T, label.size = 3, label.box = TRUE, raster = F)
p.sex
saveRDS(Seurat_object, file = "Seurat_object_raw.rds")
###########################standard#####################
# must be sequential, or will be crashed
Seurat_object <- SCTransform(Seurat_object, vars.to.regress = c("nCount_RNA", "percent.ribo", "percent.mt"), conserve.memory = F, verbose = T)
#Seurat_object <- NormalizeData(Seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
###########################highly variable genes#####################
Seurat_object <- FindVariableFeatures(Seurat_object, 
                                      selection.method = "vst", 
                                      nfeatures = 3000,
                                      mean.cutoff = c(0.003, 2), 
                                      dispersion.cutoff = c(0,1))
##########################scale#####################
Seurat_object <- ScaleData(Seurat_object, features = rownames(Seurat_object), model.use = "poisson", use.umi = TRUE)

#Seurat_object <- ScaleData(Seurat_object, vars.to.regress = c("nCount_RNA", "percent.mt"), features = rownames(Seurat_object), model.use = "poisson", use.umi = TRUE)

 #######################embeddings#####################
#Seurat_object <- RunPCA(Seurat_object, features = VariableFeatures(object = Seurat_object))
Seurat_object <- RunPCA(Seurat_object, npcs = 100)
Seurat_object <- FindNeighbors(Seurat_object, dims = 1:30)
Seurat_object <- FindClusters(Seurat_object, resolution = seq(0.1, 1.5, by=0.2))
Seurat_object <- RunUMAP(Seurat_object, dims = 1:30)
Seurat_object <- RunTSNE(Seurat_object, dims = 1:30)
saveRDS(Seurat_object, file = 'Seurat_object.rds')
###################correct batch effect##########################
library(stringr)
#batch
group_values <- ifelse(str_detect(colnames(Seurat_object), "F30|F38|F12|F14|F25"), "12F", 
                       ifelse(str_detect(colnames(Seurat_object), "M22|M17"), "1M", 
                              ifelse(str_detect(colnames(Seurat_object), "F23|F9|F18"), "2F",
                                     ifelse(str_detect(colnames(Seurat_object), "M21|M4"), "2M",
                                            ifelse(str_detect(colnames(Seurat_object), "F24|F26|F29|F32|F19|F20|F27|F28"), "3F",
                                                   ifelse(str_detect(colnames(Seurat_object), "M3|M12|M13|M6|M10|M14|M33"), "3M",
                                                          ifelse(str_detect(colnames(Seurat_object), "M9|M15|M20|M24_2|M8|M11|M30|M32"), "4M",
                                                                 ifelse(str_detect(colnames(Seurat_object), "M7|M27|M31|M1|M18|M26|M28"), "5M",
                                                                        ifelse(str_detect(colnames(Seurat_object), "M2|M16|M19|M24|M29|M5|M23|M34"), "6M",
                                                                               ifelse(str_detect(colnames(Seurat_object), "F10|F31|F33|F3|F4|F5"), "6F",
                                                                                      ifelse(str_detect(colnames(Seurat_object), "F7|F13|F21|F22|F2|F15|F16|F17"), "7F",
                                                                                             ifelse(str_detect(colnames(Seurat_object), "F34|F35|F36|F37|F1|F6|F8|F11"), "8F", NA))))))))))))
Seurat_object$batch <- group_values

## check
#b1_data <- Seurat_object[, Seurat_object$batch == "B1"]
#head(b1_data)
#dim(b1_data)
#[1] 30062  5371
Seurat_object@meta.data$batch <- as.factor(Seurat_object@meta.data$batch)
Seurat_object_harmony <- RunHarmony(Seurat_object ,"batch", plot_convergence = TRUE)
harmony_embeddings<- Embeddings(Seurat_object_harmony, "harmony")
Seurat_object_harmony <- Seurat_object_harmony %>% RunUMAP(reduction="harmony", dims = 1:30) %>% RunTSNE(reduction="harmony", dims = 1:30) %>% FindNeighbors(reduction = "harmony", dims = 1:30)
Seurat_object <- FindClusters(Seurat_object_harmony, resolution = seq(0.1, 1.5, by=0.2), algorithm = 1)
colnames(Seurat_object@meta.data)
Seurat_object@meta.data$RNA_snn_res.1.5 <- as.numeric(as.character(Seurat_object@meta.data$RNA_snn_res.1.5))
# rank
order(Seurat_object@meta.data$RNA_snn_res.1.5)
Seurat_object@meta.data$RNA_snn_res.1.5 <- factor(Seurat_object@meta.data$RNA_snn_res.1.5)
Seurat_object@meta.data$RNA_snn_res.1.5 <- factor(Seurat_object@meta.data$RNA_snn_res.1.5, levels = sort(unique(Seurat_object@meta.data$RNA_snn_res.1.5)))
levels(Seurat_object@meta.data$RNA_snn_res.1.5)
saveRDS(Seurat_object, file = 'Seurat_object.rds')

################celltypist#######################
#celltypist.ipynb
#celltypist.R
counts_matrix <- Seurat_object[["SCT"]]@counts
meta_data <- Seurat_object@meta.data
library('Matrix')
library('reticulate')
library('Seurat')
#use_condaenv(condaenv = "py3.11", required = TRUE)
use_python("/home/lfzhang/miniconda3/envs/py3.11/bin/python")
py_config()
numpy = import("numpy")
pandas = import("pandas")
mpl = import("matplotlib")
mpl$use("Agg")
scanpy = import("scanpy")
celltypist = import("celltypist")
#celltypist$models$download_models(force_update = T)
model <- celltypist$models$Model$load(model = '/home/lfzhang/.celltypist/data/models/Adult_Human_PrefrontalCortex.pkl')
print(model)
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(counts_matrix))),
                       obs = pandas$DataFrame(meta_data),
                       var = pandas$DataFrame(data.frame(gene = rownames(counts_matrix),
                                                         row.names = rownames(counts_matrix)))
)
gc()
py$adata <- adata
reticulate::py_run_string("
import numpy as np
adata.X = np.array(adata.X)
")
py_run_string("adata.X.setflags(write=True)")
reticulate::py_run_string("
import scanpy as sc
import celltypist
adata_copy = sc.AnnData(X=adata.X.copy(), obs=adata.obs, var=adata.var)
sc.pp.normalize_total(adata_copy, target_sum=1e4)
sc.pp.log1p(adata_copy)
predictions = celltypist.annotate(adata_copy, model = 'Adult_Human_PrefrontalCortex.pkl', majority_voting = True)
")
predictions <- py$predictions
Seurat_object <- readRDS("~/SingleCell/Seurat_object.rds")
Seurat_object = AddMetaData(Seurat_object, predictions$predicted_labels)
saveRDS(Seurat_object, file = 'Seurat_object.rds')

###################SingleR###################
library(SingleR)
library(celldex)
# obtain DB
hpca.se <- HumanPrimaryCellAtlasData()# create file .cache/R/ExperimentHub
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

###############doubletfinder###############
#DoubletFinder
library(DoubletFinder)
#Seurat_object <- readRDS('Seurat_object_harmony.rds')

DefaultAssay(Seurat_object) <- "SCT"
sweep.res.list <- paramSweep_v3(Seurat_object, PCs = 1:30, sct = T, num.cores = 1)
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)  
bcmvn <- find.pK(sweep.stats) # best parameter
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
# pk_bcmvn = 0.0,001

annotations <- Seurat_object@meta.data$MainType
homotypic.prop <- modelHomotypic(annotations)  
#DoubletRate = ncol(Seurat_object)*7.6*1e-6 # according to 10X Help document, doublet rates increase by 7.6% per 1000 cells
DoubletRate = 0.075 # ~

# estimate the proportion of homotypic doublets by artificially mixing doublets from seurat_clusters according to the parameters in modelHomotypic().
nExp_poi <- round(DoubletRate*length(Seurat_object$MainType)) 
# calculate doublet rate
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

Seurat_object <- doubletFinder_v3(Seurat_object, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
Seurat_object <- doubletFinder_v3(Seurat_object, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = 'pANN_0.25_0.001_10840', sct = T)

Seurat_object@meta.data[,"DF_hi.lo"] <- Seurat_object@meta.data$DF.classifications_0.25_0.001_7929
Seurat_object@meta.data$DF_hi.lo[which(Seurat_object@meta.data$DF_hi.lo == "Doublet" & Seurat_object@meta.data$DF.classifications_0.25_0.001_7929 == "Singlet")] <- "Doublet-Low Confidience"
Seurat_object@meta.data$DF_hi.lo[which(Seurat_object@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
table(Seurat_object@meta.data$DF_hi.lo)
#Doublet-High Confidience                  Singlet 
#7929                   136602 
singlet <- subset(Seurat_object, subset = DF.classifications_0.25_0.001_7929 == "Singlet")


colnames(Seurat_object@meta.data)
# rank meta data
sorted_colnames <- colnames(Seurat_object@meta.data)[order(tolower(colnames(Seurat_object@meta.data)))]
Seurat_object@meta.data <- Seurat_object@meta.data[, sorted_colnames]


## results
png("Figure/doubletFinder_new.png",2500,1800,res=300)
DimPlot(Seurat_object, reduction = "umap", group.by ="DF_hi.lo",cols =c("#1a3d85","#d76092","gold"),raster = F)
dev.off()

singlet <- subset(Seurat_object, DF.classifications_0.25_0.001_7929 == "Singlet")
saveRDS(singlet, file = "singlet.rds")
saveRDS(Seurat_object, file = "Seurat_object.rds")

Seurat_object@meta.data$doublet <- Seurat_object@meta.data$DF.classifications_0.25_0.15_8132
Idents(Seurat_object) <- Seurat_object$doublet
pdf(paste0(res_home, "Figure/Seurat_object_doubletfinder.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(2),label = T, repel = T,label.size = 3,label.box = T)
dev.off()


#DoubletDecon
library(DoubletDecon)
Idents(Seurat_object) <- Seurat_object$RNA_snn_res.0.4
newFiles <- Improved_Seurat_Pre_Process(Seurat_object, num_genes=50, write_files=FALSE)

#saveRDS(newFiles, file = "DoubletDecon_newFiles.rds")

class(as.data.frame(newFiles$newGroupsFile))


filename="Seurat_object"
location="./"
write.table(newFiles$newExpressionFile, paste0(location, filename, "_expression.txt"), sep="\t")
write.table(newFiles$newFullExpressionFile, paste0(location, filename, "_fullExpression.txt"), sep="\t")
write.table(newFiles$newGroupsFile, paste0(location, filename , "_groups.txt"), sep="\t", col.names = F)
newGroupsFile_data_frame <- as.data.frame(newFiles$newGroupsFile)

results <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                              groupsFile=as.data.frame(newFiles$newGroupsFile), 
                              filename="tmp", 
                              location="./",
                              fullDataFile=NULL, 
                              removeCC=FALSE, 
                              species="hsa", 
                              rhop=1, 
                              write=TRUE, 
                              PMF=TRUE, 
                              useFull=FALSE, 
                              heatmap=FALSE,
                              centroids=TRUE,
                              num_doubs=100, 
                              only50=FALSE,
                              min_uniq=4)
doublet_df <- as.data.frame(results$Final_doublets_groups)
singlet_df <- as.data.frame(results$Final_nondoublets_groups)
doublet_df$DoubletDecon="Doublet"
singlet_df$DoubletDecon="Singlet"
DoubletDecon_df <- rbind(doublet_df,singlet_df)
DoubletDecon_df$barcode = rownames(DoubletDecon_df)
DoubletDecon_df <- DoubletDecon_df[,c("barcode","DoubletDecon")]
rownames(DoubletDecon_df) <- DoubletDecon_df$barcode
DoubletDecon_df <- DoubletDecon_df[rownames(Seurat_object@meta.data),]
DoubletDecon_df$DoubletDecon
Seurat_object@meta.data$doublet_DoubletDecon <- DoubletDecon_df$DoubletDecon
Idents(Seurat_object) <- Seurat_object$doublet_DoubletDecon
DimPlot(Seurat_object, reduction = 'umap',cols = c("#CD9B9B", "#CD2626"),label = T, repel = T,label.size = 3,label.box = T)
singlet <- subset(Seurat_object, subset = doublet_DoubletDecon == "Singlet")
#reclustering
singlet <- NormalizeData(singlet, normalization.method = "LogNormalize", scale.factor = 10000)
singlet <- FindVariableFeatures(singlet, 
                                selection.method = "vst", 
                                nfeatures = 2000,
                                mean.cutoff = c(0.003, 2), 
                                dispersion.cutoff = c(0,1))
singlet <- ScaleData(singlet, features = rownames(singlet))
singlet <- RunPCA(singlet, features = VariableFeatures(object = singlet))
singlet <- FindNeighbors(singlet, dims = 1:80)
singlet <- FindClusters(singlet, resolution = seq(0.1, 1.5, by=0.2))
singlet <- RunUMAP(singlet, dims = 1:80)
singlet <- RunTSNE(singlet, dims = 1:80)
Seurat_object <- singlet

#############annotation################
OLG <- c("SLC44A1","PLP1","MBP", "MOG", "MOBP")
ASC <- c("SLC1A2","ALDH1L1", "ALDH1A1", "GFAP", "GLUL", "GJA1", "SOX9", "AQP4", "NDRG2")
NEURON <- c("SNAP25", "RBFOX3")
EPD <- c("CLDN5", "VIM")
OPC <- c("PTPRZ1","PCDH15","OLIG1","OLIG2","PDGFRA")
MicroM <- c("PTPRC","CSF1R","APBB1IP","P2RY12","CX3CR1","ITGAM")
Ex_1 <- c("SATB2", "SLC17A7")
Inhib <- c("GAD1", "GAD2")
sub_inhib <- c("VIP", "PVALB", "SST", "ADARB2", "LHX6", "LAMP5", "PAX6")

.cluster_cols <- c("#7d7db2","#c3556f","#f2d2cd","#e4ce84","#65a99e")
markers <- c(Ex_1,Inhib,OLG,OPC,ASC,MicroM,EPD)
cols_fun <- colorRampPalette(.cluster_cols)
cols <- cols_fun(length(table(Seurat_object$MainType)))
print(cols)
#[1] "#7D7DB2" "#A4668B" "#C9667C" "#E4AEB2" "#EED0B8" "#E6CE8E" "#ADBE8F"
#[8] "#65A99E"
.cluster_cols <- c("#7D7DB2","#C9667C","#E4AEB2","#EED0B8","#E6CE8E","#ADBE8F","#65A99E")
barplot(rep(1, length(table(Seurat_object$MainType))), col=cols, space=0)
unique(Seurat_object$MainType)
Idents(Seurat_object) <- "MainType"
.cluster_cols <- c("#7D7DB2","#C9667C","#E4AEB2","#EED0B8","#E6CE8E","#ADBE8F","#65A99E")
names(.cluster_cols) <- c("Micro","Ex","Inhib","ASC","OPC","EPD","OLG")
DimPlot(Seurat_object, group.by = "MainType",cols = .cluster_cols, pt.size = 0.5, raster = F)
Idents(Seurat_object) <- Seurat_object$SCT_snn_res.0.3
DotPlot(Seurat_object, features = markers, group.by = "SCT_snn_res.0.3") + coord_flip()+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#7D7DB2","grey","#c3556f"))
ggsave(filename = "DotplotClusters.png",path = paste0(res_home, "Figure/"),height=10, width=10)
Seurat_object@meta.data$MainType <- NA
Seurat_object$MainType[which(Seurat_object$SCT_snn_res.0.3 %in% c(2,16)) ]<- "Oli"
Seurat_object$MainType[which(Seurat_object$SCT_snn_res.0.3 %in% c(4)) ]<- "Ast"
Seurat_object$MainType[which(Seurat_object$SCT_snn_res.0.3 %in% c(6,7,9,13)) ]<- "InN"
Seurat_object$MainType[which(Seurat_object$SCT_snn_res.0.3 %in% c(1,3,5,10,11,15)) ]<- "ExN"
Seurat_object$MainType[which(Seurat_object$SCT_snn_res.0.3 %in% c(8)) ]<- "OPC"
Seurat_object$MainType[which(Seurat_object$SCT_snn_res.0.3 %in% c(12)) ]<- "Epd"
Seurat_object$MainType[which(Seurat_object$SCT_snn_res.0.3 %in% c(14)) ]<- "Mic"
Idents(Seurat_object) <- Seurat_object$MainType
pumap <- DimPlot(Seurat_object, reduction = 'umap',cols = .cluster_cols,label = T,repel = T,label.size = 4,label.box = T, raster = F)
pumap
ggsave(filename = "Seurat_object_UMAP_MainType.png", path = paste0(res_home, "Figure/"), plot = pumap,  dpi = 400, width = 6, height = 6, limitsize = F)

Seurat_object$MainType <- as.character(Seurat_object$MainType)
Seurat_object$MainType <- factor(Seurat_object$MainType, levels = c("Ex","Inhib","OLG","OPC","ASC","Micro","EPD"))

p <- DotPlot(Seurat_object, features = markers, group.by = "MainType") + coord_flip() +
  scale_color_gradientn(values = seq(0,1,0.2), colours = c("grey","#c3556f")) +
  labs(x="Markers",y="Cell Types")
p
ggsave(filename = "DotplotMainType.png",path = paste0(res_home, "Figure/"),height=8, width=6)

library(Seurat)
library(dplyr)
library(tidyr)

avg_exp <- AverageExpression(Seurat_object, group.by = "MainType")
avg_exp_df <- avg_exp$SCT
avg_exp_df <- as.data.frame(avg_exp_df)
avg_exp_df$gene <- rownames(avg_exp_df)
long_exp_df <- pivot_longer(avg_exp_df, 
                            cols = -gene, 
                            names_to = "MainType", 
                            values_to = "expression")
head(long_exp_df)

filtered_data <- long_exp_df %>%
  filter(gene %in% markers)
library(pheatmap)
heatmap_data <- filtered_data %>%
  pivot_wider(names_from = MainType, values_from = expression, values_fill = list(expression = 0))

gene_names <- heatmap_data$gene
heatmap_matrix <- as.matrix(heatmap_data[,-1])
heatmap <- pheatmap(heatmap_matrix, 
         labels_row = gene_names,
         scale = "row", 
         clustering_distance_rows = "minkowski", 
         clustering_distance_cols = "minkowski", 
         border_color = "white",
         clustering_method = "ward.D",
         color = colorRampPalette(c("#b79f53", "white", "#5d7b4e"))(100))
ggsave(heatmap, filename = "heatmap_MainType.png",path = paste0(res_home, "Figure/"),height=6, width=6)
