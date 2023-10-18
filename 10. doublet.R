#用了DoubletFinder和DoubletDecon这两个包分别测了两次看看区别
library(DoubletFinder)
Seurat_object <- readRDS('Seurat_object_harmony.rds')

sweep.res.list <- paramSweep_v3(Seurat_object, PCs = 1:17, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
# pk_bcmvn = 0.005

annotations <- Seurat_object@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(Seurat_object)*7.6*1e-6 #按每增加1000个细胞，双细胞比率增加千分之7.6来计算
DoubletRate = 0.15 # 0.4229324有点太高了，用另一个软件测出来是0.1416916746，所以定成0.15再看看


#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
nExp_poi <- round(DoubletRate*length(Seurat_object$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
# 计算双细胞比例
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

Seurat_object <- doubletFinder_v3(Seurat_object, PCs = 1:17, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
Seurat_object <- doubletFinder_v3(Seurat_object, PCs = 1:17, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)

colnames(Seurat_object@meta.data)
# 感觉元数据里的列名太混乱了，于是根据列名的首字母对列进行了排序
sorted_colnames <- colnames(Seurat_object@meta.data)[order(tolower(colnames(Seurat_object@meta.data)))]
# 使用排序后的列名重排meta.data中的列
Seurat_object@meta.data <- Seurat_object@meta.data[, sorted_colnames]



Seurat_object@meta.data[,"DF_hi.lo"] <- Seurat_object@meta.data$DF.classifications_0.25_0.26_21742
Seurat_object@meta.data$DF_hi.lo[which(Seurat_object@meta.data$DF_hi.lo == "Doublet" & Seurat_object@meta.data$DF.classifications_0.25_0.26_21742 == "Singlet")] <- "Doublet-Low Confidience"
Seurat_object@meta.data$DF_hi.lo[which(Seurat_object@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
table(Seurat_object@meta.data$DF_hi.lo)
#Doublet-High Confidience                  Singlet 
#21742                   30875 


## 结果展示，分类结果在pbmc@meta.data中
png("Figure/doubletFinder_new.png",2500,1800,res=300)
DimPlot(Seurat_object, reduction = "umap", group.by ="DF_hi.lo",cols =c("black","red","gold"))
dev.off()

saveRDS(Seurat_object, file = "Seurat_object.rds")

#选一个同源/非同源双细胞作为最终结果
Seurat_object@meta.data$doublet <- Seurat_object@meta.data$DF.classifications_0.25_0.26_21742
Idents(Seurat_object) <- Seurat_object$doublet
pdf(paste0(res_home, "Figure/Seurat_object_doubletfinder.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = colorRampPalette(.cluster_cols)(2),label = T, repel = T,label.size = 3,label.box = T)
dev.off()



library(DoubletDecon)
Idents(Seurat_object) <- Seurat_object$RNA_snn_res.0.2
newFiles <- Improved_Seurat_Pre_Process(Seurat_object, num_genes=50, write_files=FALSE)

#saveRDS(newFiles, file = "DoubletDecon_newFiles.rds")

class(as.data.frame(newFiles$newGroupsFile))


filename="Seurat_object"
location="./"
write.table(newFiles$newExpressionFile, paste0(location, filename, "_expression.txt"), sep="\t")
write.table(newFiles$newFullExpressionFile, paste0(location, filename, "_fullExpression.txt"), sep="\t")
write.table(newFiles$newGroupsFile, paste0(location, filename , "_groups.txt"), sep="\t", col.names = F)


class(newFiles$newExpressionFile)
class(newGroupsFile_data_frame)
newGroupsFile_data_frame <- as.data.frame(newFiles$newGroupsFile)
subset_expression <- newFiles$newExpressionFile[1:472, 1:100]
subset_groups <- newGroupsFile_data_frame[1:100, ]
write.table(subset_expression, "subset_expression.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(subset_groups, "subset_groups.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

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
pdf(paste0(res_home, "Figure/Seurat_object_DoubletDecon.pdf"),width = 7,height = 6)
DimPlot(Seurat_object, reduction = 'umap',cols = c("#D76092", "#008DCE"),label = T, repel = T,label.size = 3,label.box = T)
dev.off()

saveRDS(Seurat_object, file = "Seurat_object.rds")
