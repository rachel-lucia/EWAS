# redirect console output to text file
setwd("/pub/mcfarlar/")
output <- file("./GP/ComBatOutput.txt")
sink(file = output, append = TRUE, type = "output", split = TRUE)
sink(file = output, append = TRUE, type = "message")

# load previous R session
load("LoadFilter.RData")

# load required packages
#library(minfi)
library(sva)
library(ChAMP)

# setup sample sheet (substitute NA values w/ median or most frequent)
samplesRaw <- read.csv("SampleSheetGP.csv",header = TRUE, row.names = 1) 
samples <- samplesRaw
samples$Alcohol[is.na(samples$Alcohol)] <- 1
samples$Smoking[is.na(samples$Smoking)] <- 0
levels(samples$Race) <- c("Asian","Black","Hispanic","Other","Other","White")
samples$Race <- relevel(samples$Race, "White")
samples$EatOrganic[is.na(samples$EatOrganic)] <- 2
samples$HEI[is.na(samples$HEI)] <- 63.86

samples$Sample_Plate <- as.factor(samples$Sample_Plate)
samples$ColBatch <- paste(samples$Column_N,samples$Batch,sep = "_")
samples$ColBatch[samples$Batch == "Pilot"] <- "Pilot"

# setup beta matrix, removing dropped samples
myUncorrBetas <- myNormed[,colnames(myNormed) %in% rownames(samples)]
# colnames(myUncorrBetas) == rownames(samples) # yup

# convert to M values
myMvals <- log2(myUncorrBetas/(1-myUncorrBetas))

# visualize batch effects on M values
samples$logGP = log(samples$GP)
samples$logAMPA = log(samples$AMPA)
keep <- c("Sentrix_ID","Batch","Sample_Plate","Column_N","Row",
          "Age","Race","Alcohol","Smoking","BMI",
          "EatOrganic","HEI","logGP","logAMPA")
SVDdata <- samples[,keep]
champ.SVD(beta=myMvals, pd=SVDdata, resultsDir = "./SVD/") # takes ~30m

# run ComBat to correct batch effects, protecting clinical variables
## first var to adjust for: Sentrix_ID (accounts for row and plate)
batch = as.factor(samples$Sentrix_ID)
# this round, full model with var of interest and all covariates
myModel <- model.matrix(~ Race + BMI + Age + Alcohol + Smoking + EatOrganic + HEI + logGP + logAMPA, data = samples)
myAdjMs = ComBat(dat=myMvals, batch=batch, mod=myModel) # 3 batches = 6 mins; 51 = 30 mins

## do other batch effects persist? if so, run ComBat again
# visualize - run Champ.SVD again in M-space
champ.SVD(beta=myAdjMs, pd=SVDdata, resultsDir = "./SVD/Post1/") # takes ~30m

# column is still in top 5 batches + still an effect of batch. adjust for this as well
batch2 = as.factor(samples$ColBatch)
myAdjMs2 = ComBat(dat=myAdjMs, batch=batch2, mod=myModel)

# final SVD
champ.SVD(beta=myAdjMs2, pd=SVDdata, resultsDir = "./SVD/Post2/") # takes ~30m


# check on this
rm(myAdjBs,myMvals,myFilter,myNormed,myAdjMs,myUncorrBetas)
gc()

# save R session, commands, and filtered/normalized beta values for later
sessionInfo() # saved below, ran this interactively
save(file="./GP/ComBatOut.RData",list = ls(all.names = TRUE), envir = .GlobalEnv)
q(save = "no", status = 0)


# > sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS release 6.9 (Final)
# 
# Matrix products: default
# BLAS: /usr/lib64/libblas.so.3.2.1
# LAPACK: /usr/lib64/atlas/liblapack.so.3.0
# 
# locale:
#   [1] LC_CTYPE=en_US.iso885915       LC_NUMERIC=C
# [3] LC_TIME=en_US.iso885915        LC_COLLATE=en_US.iso885915
# [5] LC_MONETARY=en_US.iso885915    LC_MESSAGES=en_US.iso885915
# [7] LC_PAPER=en_US.iso885915       LC_NAME=C
# [9] LC_ADDRESS=C                   LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.iso885915 LC_IDENTIFICATION=C
# 
# attached base packages:
#   [1] splines   stats4    parallel  stats     graphics  grDevices utils
# [8] datasets  methods   base
# 
# other attached packages:
#   [1] ChAMP_2.12.4
# [2] IlluminaHumanMethylationEPICmanifest_0.3.0
# [3] Illumina450ProbeVariants.db_1.18.0
# [4] DMRcate_1.18.0
# [5] DMRcatedata_1.18.0
# [6] DSS_2.30.1
# [7] bsseq_1.18.0
# [8] FEM_3.10.0
# [9] graph_1.60.0
# [10] org.Hs.eg.db_3.7.0
# [11] impute_1.56.0
# [12] igraph_1.2.4.1
# [13] corrplot_0.84
# [14] marray_1.60.0
# [15] limma_3.38.3
# [16] Matrix_1.2-17
# [17] AnnotationDbi_1.44.0
# [18] ChAMPdata_2.14.1
# [19] minfi_1.28.4
# [20] bumphunter_1.24.5
# [21] locfit_1.5-9.1
# [22] iterators_1.0.10
# [23] foreach_1.4.4
# [24] Biostrings_2.50.2
# [25] XVector_0.22.0
# [26] SummarizedExperiment_1.12.0
# [27] DelayedArray_0.8.0
# [28] matrixStats_0.54.0
# [29] Biobase_2.42.0
# [30] GenomicRanges_1.34.0
# [31] GenomeInfoDb_1.18.2
# [32] IRanges_2.16.0
# [33] S4Vectors_0.20.1
# [34] BiocGenerics_0.28.0
# [35] sva_3.30.1
# [36] BiocParallel_1.16.6
# [37] genefilter_1.64.0
# [38] mgcv_1.8-28
# [39] nlme_3.1-140
# 
# loaded via a namespace (and not attached):
#   [1] R.utils_2.8.0
# [2] tidyselect_0.2.5
# [3] IlluminaHumanMethylationEPICanno.ilm10b2.hg19_0.6.0
# [4] RSQLite_2.1.1
# [5] htmlwidgets_1.3
# [6] combinat_0.0-8
# [7] grid_3.5.1
# [8] munsell_0.5.0
# [9] codetools_0.2-16
# [10] preprocessCore_1.44.0
# [11] nleqslv_3.3.2
# [12] statmod_1.4.30
# [13] withr_2.1.2
# [14] fastICA_1.2-1
# [15] colorspace_1.4-1
# [16] knitr_1.23
# [17] rstudioapi_0.10
# [18] kpmt_0.1.0
# [19] JADE_2.0-1
# [20] isva_1.9
# [21] GenomeInfoDbData_1.2.0
# [22] bit64_0.9-7
# [23] rhdf5_2.26.2
# [24] xfun_0.7
# [25] biovizBase_1.30.1
# [26] doParallel_1.0.14
# [27] R6_2.4.0
# [28] clue_0.3-57
# [29] illuminaio_0.24.0
# [30] AnnotationFilter_1.6.0
# [31] bitops_1.0-6
# [32] reshape_0.8.8
# [33] assertthat_0.2.1
# [34] promises_1.0.1
# [35] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0
# [36] scales_1.0.0
# [37] nnet_7.3-12
# [38] gtable_0.3.0
# [39] affy_1.60.0
# [40] methylumi_2.28.0
# [41] ensembldb_2.6.8
# [42] rlang_0.3.4
# [43] rtracklayer_1.42.2
# [44] lazyeval_0.2.2
# [45] acepack_1.4.1
# [46] GEOquery_2.50.5
# [47] dichromat_2.0-0
# [48] checkmate_1.9.3
# [49] reshape2_1.4.3
# [50] BiocManager_1.30.4
# [51] GenomicFeatures_1.34.8
# [52] backports_1.1.4
# [53] httpuv_1.5.1
# [54] qvalue_2.14.1
# [55] Hmisc_4.2-0
# [56] tools_3.5.1
# [57] nor1mix_1.2-3
# [58] ggplot2_3.1.1
# [59] affyio_1.52.0
# [60] lumi_2.34.0
# [61] RColorBrewer_1.1-2
# [62] DNAcopy_1.56.0
# [63] siggenes_1.56.0
# [64] Rcpp_1.0.1
# [65] plyr_1.8.4
# [66] base64enc_0.1-3
# [67] ROC_1.58.0
# [68] progress_1.2.2
# [69] zlibbioc_1.28.0
# [70] purrr_0.3.2
# [71] RCurl_1.95-4.12
# [72] BiasedUrn_1.07
# [73] prettyunits_1.0.2
# [74] rpart_4.1-15
# [75] openssl_1.3
# [76] viridis_0.5.1
# [77] cluster_2.0.9
# [78] magrittr_1.5
# [79] data.table_1.12.2
# [80] goseq_1.34.1
# [81] wateRmelon_1.26.0
# [82] ProtGenerics_1.14.0
# [83] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
# [84] missMethyl_1.16.0
# [85] RPMM_1.25
# [86] evaluate_0.13
# [87] mime_0.6
# [88] hms_0.4.2
# [89] xtable_1.8-4
# [90] globaltest_5.36.0
# [91] XML_3.98-1.19
# [92] mclust_5.4.3
# [93] gridExtra_2.3
# [94] compiler_3.5.1
# [95] biomaRt_2.38.0
# [96] tibble_2.1.1
# [97] KernSmooth_2.23-15
# [98] crayon_1.3.4
# [99] R.oo_1.22.0
# [100] htmltools_0.3.6
# [101] later_0.8.0
# [102] Formula_1.2-3
# [103] tidyr_0.8.3
# [104] DBI_1.0.0
# [105] geneLenDataBase_1.18.0
# [106] MASS_7.3-51.4
# [107] readr_1.3.1
# [108] permute_0.9-5
# [109] quadprog_1.5-7
# [110] R.methodsS3_1.7.1
# [111] Gviz_1.26.5
# [112] pkgconfig_2.0.2
# [113] GenomicAlignments_1.18.1
# [114] registry_0.5-1
# [115] IlluminaHumanMethylation450kmanifest_0.4.0
# [116] foreign_0.8-71
# [117] plotly_4.9.0
# [118] xml2_1.2.0
# [119] annotate_1.60.1
# [120] rngtools_1.3.1
# [121] pkgmaker_0.27
# [122] multtest_2.38.0
# [123] beanplot_1.2
# [124] ruv_0.9.7
# [125] bibtex_0.4.2
# [126] doRNG_1.7.1
# [127] stringr_1.4.0
# [128] VariantAnnotation_1.28.13
# [129] digest_0.6.19
# [130] rmarkdown_1.13
# [131] base64_2.0
# [132] htmlTable_1.13.1
# [133] dendextend_1.12.0
# [134] DelayedMatrixStats_1.4.0
# [135] curl_3.3
# [136] shiny_1.3.2
# [137] Rsamtools_1.34.1
# [138] gtools_3.8.1
# [139] jsonlite_1.6
# [140] Rhdf5lib_1.4.3
# [141] viridisLite_0.3.0
# [142] askpass_1.1
# [143] BSgenome_1.50.0
# [144] pillar_1.4.0
# [145] lattice_0.20-38
# [146] httr_1.4.0
# [147] survival_2.44-1.1
# [148] GO.db_3.7.0
# [149] glue_1.3.1
# [150] shinythemes_1.1.2
# [151] bit_1.1-14
# [152] stringi_1.4.3
# [153] HDF5Array_1.10.1
# [154] blob_1.1.1
# [155] latticeExtra_0.6-28
# [156] memoise_1.1.0
# [157] dplyr_0.8.1
