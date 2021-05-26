# manuscript-v1

## Overview of the single cell RNA-seq analysis

Experiment was done in three batches:
Batch1 samples: NM05N, NM04N1, NH64T
Batch2 samples: NH87ND, NH87T, HCM-CSHL-0655-C50
Batch3 samples: multiplexed 4 samples with the following barcodes
  GACAGTGC HCM-CSHL-0366-C50
  GAGTTAGC NH85TSc
  GATGAATC NH95T
  GCCAAGAC NH93T

10x Single cell RNA-seq analysis
-- pre-processing using cell ranger
-- aligned to grch38

sources and tutorials used for analysis:
https://github.com/hbctraining/scRNA-seq/blob/master/lessons/03_SC_quality_control-setup.md 
https://satijalab.org/seurat/v3.0/merge_vignette.html
https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html + seurat vignette
  cluster specific marker expression from Seurat's vignette.
 
  
### For Figure 8
Following are the detail for each code: 
integrated_tumor_only_processing.Rmd: processing of individual TNBC tumor h5 files for integration. Removed low quality cells and genes

integrated_tumor_only_SCT_integration.Rmd: SC transformation followed by anchor dependent integration for the TNBC tumors



### R session information pertaining to this analysis

R version 3.6.3 (2020-02-29)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggrastr_0.2.3      plotly_4.9.3       reshape2_1.4.4     dittoSeq_1.3.6     sctransform_0.3.2  metap_1.4          RCurl_1.98-1.2    
 [8] forcats_0.5.1      purrr_0.3.4        tidyr_1.1.3        tibble_3.1.0       tidyverse_1.3.0    readr_1.4.0        ggplot2_3.3.3     
[15] stringr_1.4.0      hdf5r_1.3.3        cowplot_1.1.1      future.apply_1.7.0 future_1.21.0      magrittr_2.0.1     Matrix_1.3-2      
[22] dplyr_1.0.5        Seurat_3.2.3       knitr_1.31        

loaded via a namespace (and not attached):
  [1] readxl_1.3.1                backports_1.2.1             sn_1.6-2                    plyr_1.8.6                 
  [5] igraph_1.2.6                lazyeval_0.2.2              splines_3.6.3               BiocParallel_1.20.1        
  [9] listenv_0.8.0               scattermore_0.7             GenomeInfoDb_1.22.1         TH.data_1.0-10             
 [13] digest_0.6.27               htmltools_0.5.1.1           fansi_0.4.2                 tensor_1.5                 
 [17] cluster_2.1.1               ROCR_1.0-11                 globals_0.14.0              modelr_0.1.8               
 [21] matrixStats_0.58.0          sandwich_3.0-0              colorspace_2.0-0            rvest_1.0.0                
 [25] ggrepel_0.9.1               haven_2.3.1                 rbibutils_2.0               xfun_0.21                  
 [29] crayon_1.4.1                jsonlite_1.7.2              spatstat_1.64-1             spatstat.data_2.0-0        
 [33] survival_3.2-7              zoo_1.8-9                   glue_1.4.2                  polyclip_1.10-0            
 [37] gtable_0.3.0                zlibbioc_1.32.0             XVector_0.26.0              leiden_0.3.7               
 [41] DelayedArray_0.12.3         SingleCellExperiment_1.8.0  BiocGenerics_0.32.0         abind_1.4-5                
 [45] scales_1.1.1                pheatmap_1.0.12             mvtnorm_1.1-1               DBI_1.1.1                  
 [49] miniUI_0.1.1.1              Rcpp_1.0.6                  plotrix_3.8-1               viridisLite_0.3.0          
 [53] xtable_1.8-4                tmvnsim_1.0-2               reticulate_1.18             bit_4.0.4                  
 [57] rsvd_1.0.3                  stats4_3.6.3                htmlwidgets_1.5.3           httr_1.4.2                 
 [61] RColorBrewer_1.1-2          TFisher_0.2.0               ellipsis_0.3.1              ica_1.0-2                  
 [65] farver_2.1.0                pkgconfig_2.0.3             uwot_0.1.10                 dbplyr_2.1.0               
 [69] deldir_0.2-10               utf8_1.1.4                  labeling_0.4.2              tidyselect_1.1.0           
 [73] rlang_0.4.10                later_1.1.0.1               munsell_0.5.0               cellranger_1.1.0           
 [77] tools_3.6.3                 cli_2.3.1                   generics_0.1.0              mathjaxr_1.4-0             
 [81] broom_0.7.5                 ggridges_0.5.3              evaluate_0.14               fastmap_1.1.0              
 [85] yaml_2.2.1                  goftest_1.2-2               bit64_4.0.5                 fs_1.5.0                   
 [89] fitdistrplus_1.1-3          RANN_2.6.1                  pbapply_1.4-3               nlme_3.1-152               
 [93] mime_0.10                   xml2_1.3.2                  compiler_3.6.3              rstudioapi_0.13            
 [97] beeswarm_0.3.1              png_0.1-7                   spatstat.utils_2.0-0        reprex_1.0.0               
[101] stringi_1.5.3               RSpectra_0.16-0             lattice_0.20-41             multtest_2.42.0            
[105] vctrs_0.3.6                 mutoss_0.1-12               pillar_1.5.1                lifecycle_1.0.0            
[109] Rdpack_2.1.1                lmtest_0.9-38               RcppAnnoy_0.0.18            data.table_1.14.0          
[113] bitops_1.0-6                irlba_2.3.3                 GenomicRanges_1.38.0        httpuv_1.5.5               
[117] patchwork_1.1.1             R6_2.5.0                    promises_1.2.0.1            KernSmooth_2.23-18         
[121] gridExtra_2.3               vipor_0.4.5                 IRanges_2.20.2              parallelly_1.23.0          
[125] codetools_0.2-18            MASS_7.3-53.1               assertthat_0.2.1            SummarizedExperiment_1.16.1
[129] withr_2.4.1                 mnormt_2.0.2                GenomeInfoDbData_1.2.2      S4Vectors_0.24.4           
[133] multcomp_1.4-16             mgcv_1.8-34                 parallel_3.6.3              hms_1.0.0                  
[137] grid_3.6.3                  rpart_4.1-15                rmarkdown_2.7               Rtsne_0.15                 
[141] Biobase_2.46.0              numDeriv_2016.8-1.1         shiny_1.6.0                 lubridate_1.7.10           
[145] ggbeeswarm_0.6.0           
