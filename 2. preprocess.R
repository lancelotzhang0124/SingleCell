setwd("/")#你的工作路径
##########调用多线程(可选)#####################
library("future")
library("future.apply")
plan(multisession)
#设置最大分配内存大小，我这里是128G
options(future.globals.maxSize = 128 * 1024^3)
future::plan(future::sequential)
#########构建表达矩阵，核心就是得到一个基因×细胞的矩阵#####################
matrix_dir = "/"
#根据实际文件夹进行修改，我这因为数据集没直接给表达矩阵，使用的是本地导入，可以直接用GEOquery包导入
barcode.path <- paste0(matrix_dir, "CellNames.csv.gz")
features.path <- paste0(matrix_dir, "GeneNames.csv.gz")
matrix.path <- paste0(matrix_dir, "GeneBarcodeMatrix_Annotated.mtx")
mat <- readMM(file = matrix.path)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names <- barcode.names[-1, ]#我这第一列是有序号的，所以都去掉了
barcode.names <- sub("^\\d+,", "", barcode.names)#这里是把分隔的内容去掉了好像，不记得了
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
feature.names <- feature.names[-1, ]
feature.names <- sub("^\\d+,", "", feature.names)
colnames(mat) = barcode.names
rownames(mat) = feature.names
#以上都是一些预处理，可以用txt文档打开先看看数据表格是什么格式的，然后按对应的方式处理
