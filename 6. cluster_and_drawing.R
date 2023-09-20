#########################降维聚类#####################
Seurat_object <- RunPCA(Seurat_object, features = VariableFeatures(object = Seurat_object))
Seurat_object <- FindNeighbors(Seurat_object, dims = 1:30)
Seurat_object <- FindClusters(Seurat_object, resolution = 0.5)
Seurat_object <- RunUMAP(Seurat_object, dims = 1:30)
Seurat_object <- RunTSNE(Seurat_object, dims = 1:30)

#########################cluster作图#####################
# 根据暴露提取组别
group_values <- ifelse(str_detect(colnames(Seurat_object), "Control"), "Control", 
                       ifelse(str_detect(colnames(Seurat_object), "Suicide"), "Suicide", NA))
group_values
Seurat_object$group <- group_values


# 使用新创建的'group'列绘制图形
p1 <- DimPlot(Seurat_object, reduction = "tsne", group.by = "group", 
              pt.size=0.8) + 
  theme(axis.text.x = element_text(size = 30),   # x轴的字体大小
        axis.text.y = element_text(size = 30), 
        legend.text = element_text(size = 30),
        legend.position = c(0.85, 0.1)
  )
p1 <- p1 + labs(title = "Group") + 
  theme(plot.title = element_text(size = 50))
ggsave("p1.pdf", p1, width=20, height=20, limitsize = FALSE)

p2 <- DimPlot(Seurat_object, reduction = "tsne", group.by = "seurat_clusters", 
              pt.size=0.8, label = TRUE, repel = TRUE, label.size = 10) + 
  theme(axis.text.x = element_text(size = 25),   # x轴的字体大小
        axis.text.y = element_text(size = 25), 
        legend.text = element_text(size = 25),
        legend.position = "right"
  ) 
p2 <- p2 + 
  coord_fixed(ratio = 1) + 
  labs(title = "Cell Type") + 
  theme(
    plot.title = element_text(size = 50, hjust = 0.5),
    plot.title.position = "plot"
  ) + 
  guides(color = guide_legend(override.aes = list(size = 6), ncol = 2))
ggsave("p2.pdf", p2, width=20, height=20, limitsize = FALSE)

p3 <- DimPlot(Seurat_object, reduction = "umap", group.by = "group", 
              pt.size=0.8)+ 
  theme(axis.text.x = element_text(size = 30),   # x轴的字体大小
        axis.text.y = element_text(size = 30), 
        legend.text = element_text(size = 30),
        legend.position = c(0.85, 0.1)
  )
p3 <- p3 + labs(title = "Group") + 
  theme(plot.title = element_text(size = 50))
ggsave("p3.pdf", p3, width=20, height=20, limitsize = FALSE)

p4 <- DimPlot(Seurat_object, reduction = "umap", group.by = "seurat_clusters", 
              pt.size=0.8, label = TRUE, repel = TRUE, label.size = 8) + 
  theme(axis.text.x = element_text(size = 30),   # x轴的字体大小
        axis.text.y = element_text(size = 30), 
        legend.text = element_text(size = 30),
        legend.position = "right"
  )
p4 <- p4 + 
  coord_fixed(ratio = 1) + 
  labs(title = "Cell Type") + 
  theme(
    plot.title = element_text(size = 50, hjust = 0.5),
    plot.title.position = "plot"
  ) + 
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggsave("p4.pdf", p4, width=20, height=20, limitsize = FALSE)
#确定维度
#Seurat_object <- JackStraw(Seurat_object, num.replicate = 100)
#Seurat_object <- ScoreJackStraw(Seurat_object, dims = 1:20)
#JackStrawPlot(Seurat_object, dims = 1)
ElbowPlot(Seurat_object)#JackStraw数据量太大的时候不适合使用
#根据这一步确定PCAdims=1:?，我确定的是到15
#PCA
VizDimLoadings(Seurat_object, dims = 1:15, reduction = "pca")
p5 <- DimPlot(Seurat_object, reduction = "pca")
ggsave("p5.pdf", p5, width=10, height=10, limitsize = FALSE )
#################为cluster进行注释#####################
all.markers <- FindAllMarkers(Seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#用鉴定出的基因替换掉我写的ur_Gene_name_0、ur_Gene_name_1
#VlnPlot(Seurat_object, features = c("ur_Gene_name_0", "ur_Gene_name_1"), slot = "counts", log = TRUE)

pdf(file="0.markerViolin01.pdf",width=10,height=6)
VlnPlot(object = Seurat_object, features = c("ur_Gene_name_0", "ur_Gene_name_1"))
dev.off()

pdf(file="0.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = Seurat_object, features = c("ur_Gene_name_0", "Aur_Gene_name_1"))
dev.off()

FeaturePlot(Seurat_object, features = c("ur_Gene_name_0"))

pdf(file="0.markerBubble.pdf",width=12,height=6)
#cluster0Marker=c("Cbr2", "Lyve1", "Selenop", "Folr2", "Ednrb", "F13a1", "Mrc1", "Igf1", "Slc40a1", "Cd163")#此处是示例，请自行更改
DotPlot(object = Seurat_object, features = cluster0Marker)
dev.off()

library("multtest")
library("metap")
DefaultAssay(Seurat_object) <- "RNA"
cluster1.markers <- FindMarkers(Seurat_object, ident.1 = 1, min.pct = 0.25)
#每个聚类前10个差异基因表达热图(如果小于10，则绘制所有标记)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
heatmap_plot <- DoHeatmap(Seurat_object, features = top10$gene) + NoLegend()
ggsave(filename = "heatmap_plot.png", plot = heatmap_plot, width = 20, height = 20)
