## ----eval=FALSE---------------------------------------------------------------
#  library("recount3")
#  library("dplyr")
#  library("stringr")
#  library("SummarizedExperiment")
#  library("tidyr")
#  library("tibble")
#  library("magrittr")

## ----eval=FALSE---------------------------------------------------------------
#  
#  # download data recount3
#  proj_info <- available_projects() %>%
#    filter(project == "SRP093386")
#  se <- create_rse(proj_info)
#  
#  count_mat <- assay(se)
#  rownames(count_mat) <- str_replace(rownames(count_mat), "\\.[:number:]+$", "")
#  
#  meta <- colData(se) %>%
#    as_tibble() %>%
#    separate(sra.sample_title, into = c("cell_line", "treatment", "mutation", "replicate"), sep = "-") %>%
#    select(cell_line, treatment, mutation, replicate)
#  
#  colnames(count_mat) <- paste0(meta$treatment, "_", meta$mutation, "_", meta$replicate)
#  
#  # subset to cell line T47D
#  count_mat <- count_mat[, meta$cell_line == "T47D"]
#  meta <- meta[meta$cell_line == "T47D", c("treatment", "mutation", "replicate")]
#  
#  T47D <- make_dds(count_mat, meta, ah_record = "AH89426")
#  T47D <- T47D[rowSums(assay(T47D))>0,]

## ----eval=FALSE---------------------------------------------------------------
#  dds <- T47D
#  dds <- filter_genes(dds, min_count = 5, min_rep = 4)
#  dds$mutation <- as.factor(dds$mutation)
#  dds$treatment <- as.factor(dds$treatment)
#  design(dds) <- ~ mutation + treatment
#  
#  # to not run DESeq2 in the main vignette,
#  # wo pre-compute the dispersion plot and diff testing results
#  dds <- DESeq(dds, parallel=T)
#  
#  png(filename="disp_ests.png", width=7, height=5, units="in", res=200)
#  plotDispEsts(dds)
#  dev.off()
#  
#  T47D_diff_testing <- lfcShrink(dds, coef = "mutation_WT_vs_D538G", lfcThreshold = log2(1.5), type = "normal", parallel = TRUE)
#  T47D_diff_testing$stat <- NULL
#  T47D_diff_testing$lfcSE <- NULL
#  T47D_diff_testing$pvalue <- NULL

