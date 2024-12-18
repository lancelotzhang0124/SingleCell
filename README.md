# SingleCell RNA-seq Pipeline

Welcome to the **SingleCell RNA-seq Pipeline** repository. This pipeline is designed for processing and analyzing single-cell RNA sequencing (scRNA-seq) data efficiently and effectively.

## Environment Requirements

- **Platform**: R Studio Server / Visual Studio Code
- **R Version**: [4.3.3](https://cran.rstudio.com/bin/windows/base/old/4.3.3/)
- **Bioconductor Version**: [3.17](https://bioconductor.org/news/bioc_3_17_release/)
- **Python Version**: [3.11](https://www.python.org/downloads/release/python-3110/)

## Key Packages

- **[Seurat](https://satijalab.org/seurat/)**: Version 4.4.0, [Hao*, Hao*, et al., Cell 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)00583-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421005833%3Fshowall%3Dtrue)
- **[scDRS](https://github.com/martinjzhang/scDRS)**: Version 1.0.2, [Zhang*, Hou*, et al. "Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data"](https://www.nature.com/articles/s41588-022-01167-z), Nature Genetics, 2022.
- **[DoubletFinder](https://github.com/ekernf01/DoubletFinder)**: Version 2.0.3, [McGinnis, C. S., et al. "DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors"](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0),  Cell Systems, 2019.
- **[geneNMF](https://github.com/carmonalab/GeneNMF)**: Version 0.6.0, [bioRxiv 2024 10.1101/2024.05.31.596823](https://www.biorxiv.org/content/10.1101/2024.05.31.596823v1)
- **[slingshot](https://github.com/kstreet13/slingshot)**: Version 2.10.0, [Street, K et al., "Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics."](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4772-0), BMC Genomics, 2018.
- **[CellChat](https://github.com/jinworks/CellChat)**: Version 2.1.0, [Jin et al. "CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics"](https://www.nature.com/articles/s41596-024-01045-4), Nature Protocols, 2024.
- Additional R and Python packages will be automatically installed as dependencies during setup.

## Contact

For questions or support, feel free to reach out via email: [zhanglingfeng@whu.edu.cn](mailto:zhanglingfeng@whu.edu.cn).

## Details
```R
> sessionInfo()
R version 4.3.1 (2023-06-16)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Ubuntu 18.04.6 LTS

Matrix products: default
BLAS/LAPACK: /home/lfZhang/miniconda3/envs/R4.3/lib/libopenblasp-r0.3.24.so;  LAPACK version 3.11.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=zh_CN.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C       

time zone: Asia/Shanghai
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] org.Hs.eg.db_3.18.0   AnnotationDbi_1.64.1  IRanges_2.36.0       
 [4] S4Vectors_0.40.2      Biobase_2.62.0        BiocGenerics_0.48.1  
 [7] clusterProfiler_4.8.3 fgsea_1.28.0          msigdbr_7.5.1        
[10] GeneNMF_0.6.2         SeuratObject_4.1.3    Seurat_4.3.0.1       

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22        splines_4.3.1           later_1.3.2            
  [4] ggplotify_0.1.2         bitops_1.0-7            tibble_3.2.1           
  [7] polyclip_1.10-7         lifecycle_1.0.4         globals_0.16.3         
 [10] lattice_0.21-8          MASS_7.3-60             dendextend_1.19.0      
 [13] SnowballC_0.7.1         magrittr_2.0.3          plotly_4.10.4          
 [16] httpuv_1.6.15           sctransform_0.4.1       sp_2.1-4               
 [19] spatstat.sparse_3.1-0   reticulate_1.38.0       cowplot_1.1.3          
 [22] pbapply_1.7-2           DBI_1.1.3               RColorBrewer_1.1-3     
 [25] abind_1.4-8             zlibbioc_1.48.0         Rtsne_0.17             
 [28] purrr_1.0.2             ggraph_2.1.0            RCurl_1.98-1.13        
 [31] yulab.utils_0.1.5       tweenr_2.0.3            GenomeInfoDbData_1.2.10
 [34] enrichplot_1.22.0       ggrepel_0.9.6           irlba_2.3.5.1          
 [37] listenv_0.9.1           spatstat.utils_3.0-5    tidytree_0.4.6         
 [40] pheatmap_1.0.12         goftest_1.2-3           spatstat.random_3.2-3  
 [43] fitdistrplus_1.1-11     parallelly_1.38.0       leiden_0.4.3.1         
 [46] codetools_0.2-19        ggforce_0.4.2           DOSE_3.28.1            
 [49] tidyselect_1.2.1        aplot_0.2.3             farver_2.1.2           
 [52] viridis_0.6.5           matrixStats_1.3.0       spatstat.explore_3.2-7 
 [55] jsonlite_1.8.8          tidygraph_1.3.0         progressr_0.14.0       
 [58] ggridges_0.5.6          survival_3.5-7          tools_4.3.1            
 [61] treeio_1.26.0           ica_1.0-3               Rcpp_1.0.13-1          
 [64] glue_1.8.0              gridExtra_2.3           qvalue_2.34.0          
 [67] GenomeInfoDb_1.38.1     dplyr_1.1.4             withr_3.0.2            
 [70] fastmap_1.2.0           fansi_1.0.6             RcppML_0.3.7           
 [73] digest_0.6.37           gridGraphics_0.5-1      R6_2.5.1               
 [76] mime_0.12               colorspace_2.1-1        scattermore_1.2        
 [79] GO.db_3.18.0            tensor_1.5              spatstat.data_3.1-2    
 [82] RSQLite_2.3.1           utf8_1.2.4              tidyr_1.3.1            
 [85] generics_0.1.3          data.table_1.16.2       graphlayouts_1.1.0     
 [88] httr_1.4.7              htmlwidgets_1.6.4       scatterpie_0.2.3       
 [91] uwot_0.2.2              pkgconfig_2.0.3         gtable_0.3.6           
 [94] blob_1.2.4              lmtest_0.9-40           XVector_0.42.0         
 [97] shadowtext_0.1.3        htmltools_0.5.8.1       scales_1.3.0           
[100] png_0.1-8               ggfun_0.1.5             reshape2_1.4.4         
[103] nlme_3.1-163            zoo_1.8-12              cachem_1.1.0           
[106] stringr_1.5.1           KernSmooth_2.23-22      parallel_4.3.1         
[109] miniUI_0.1.1.1          HDO.db_0.99.1           pillar_1.9.0           
[112] grid_4.3.1              vctrs_0.6.5             RANN_2.6.1             
[115] lsa_0.73.3              promises_1.3.0          xtable_1.8-4           
[118] cluster_2.1.4           cli_3.6.3               compiler_4.3.1         
[121] rlang_1.1.4             crayon_1.5.3            future.apply_1.11.2    
[124] labeling_0.4.3          plyr_1.8.9              fs_1.6.4               
[127] stringi_1.8.4           viridisLite_0.4.2       deldir_2.0-4           
[130] BiocParallel_1.36.0     babelgene_22.9          munsell_0.5.1          
[133] Biostrings_2.70.1       lazyeval_0.2.2          spatstat.geom_3.2-9    
[136] GOSemSim_2.28.0         Matrix_1.6-1.1          patchwork_1.3.0        
[139] bit64_4.0.5             future_1.34.0           ggplot2_3.5.0.9000     
[142] KEGGREST_1.42.0         shiny_1.9.1             ROCR_1.0-11            
[145] igraph_2.0.3            memoise_2.0.1           ggtree_3.10.0          
[148] fastmatch_1.1-4         bit_4.0.5               downloader_0.4         
[151] gson_0.1.0              ape_5.8     
```
