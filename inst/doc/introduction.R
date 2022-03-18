## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  max.print = 100,
  cache = FALSE,
  fig.width = 7,
  fig.height = 5,
  dev="png"
)

## ---- results="hide"----------------------------------------------------------
library("RNAseqQC")
library("DESeq2")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("tibble")
library("magrittr")

## -----------------------------------------------------------------------------
count_mat <- counts(T47D)
meta <- data.frame(colData(T47D))

# count matrix; rownames must be ENSEMBL gene IDs
count_mat[head(which(rowSums(count_mat) > 0)), 1:10]

# metadata of the samples, where row i corresponds to column i in the count matrix
meta

## ---- eval=FALSE--------------------------------------------------------------
#  dds <- make_dds(counts = count_mat, metadata = meta, ah_record = "AH89426")

## ---- include=FALSE-----------------------------------------------------------
dds <- T47D

## ---- eval=FALSE--------------------------------------------------------------
#  mcols(AnnotationHub()) %>%
#    as_tibble(rownames = "record_id") %>%
#    dplyr::filter(rdataclass == "EnsDb") %>%
#    dplyr::filter(str_detect(title, "sapiens"))

## -----------------------------------------------------------------------------
plot_total_counts(dds)

## -----------------------------------------------------------------------------
plot_library_complexity(dds)

## -----------------------------------------------------------------------------
plot_gene_detection(dds)

## -----------------------------------------------------------------------------
plot_biotypes(dds)

## -----------------------------------------------------------------------------
dds <- filter_genes(dds, min_count = 5, min_rep = 4)

## -----------------------------------------------------------------------------
vsd <- vst(dds)
mean_sd_plot(vsd)

## -----------------------------------------------------------------------------
map(c("1", "5", "14", "X"), ~plot_chromosome(vsd, .x))

## ---- fig.width=8, fig.height=12, out.width="95%"-----------------------------
# define new grouping variable
colData(vsd)$trt_mut <- paste0(colData(vsd)$treatment, "_", colData(vsd)$mutation)

ma_plots <- plot_sample_MAs(vsd, group = "trt_mut")
cowplot::plot_grid(plotlist = ma_plots[17:24], ncol = 2)

## -----------------------------------------------------------------------------
# set seed to control random annotation colors
set.seed(1)
plot_sample_clustering(vsd, anno_vars = c("treatment", "mutation", "replicate"), distance = "euclidean")

## -----------------------------------------------------------------------------
plot_pca(vsd, PC_x = 1, PC_y = 2, color_by = "treatment", shape_by = "mutation")

## -----------------------------------------------------------------------------
pca_res <- plot_pca(vsd, show_plot = FALSE)
plot_loadings(pca_res, PC = 1, annotate_top_n = 5)
plot_loadings(pca_res, PC = 1, highlight_genes = c("CD34", "FLT1", "MAPT"))
plot_loadings(pca_res, PC = 4, color_by = "gene_biotype", show_plot = F)$plot +
  theme(legend.position = "bottom")

## -----------------------------------------------------------------------------
dds <- estimateSizeFactors(dds)
plot_gene("CLEC2B", dds, x_var = "mutation", color_by = "treatment")
# modify the plot
plot_gene("CLEC2B", dds, x_var = "mutation", color_by = "treatment", show_plot = F)$plot + ggsci::scale_color_jco()
# a custom plot type
plot_data <- plot_gene("CLEC2B", dds, show_plot = F)$data
plot_data %>%
  ggplot(aes(treatment, count, color = mutation)) +
  geom_point() +
  geom_path(
    aes(x = treatment, y = mean_count, color = mutation, group = mutation),
    data = plot_data %>%
      group_by(treatment, mutation) %>%
      summarize(mean_count = mean(count), .groups = "drop")
  ) +
  labs(y = "log2(norm. count)", title="CLEC2B") +
  cowplot::theme_cowplot()

## ---- eval=FALSE--------------------------------------------------------------
#  plots <- rownames(dds)[1:100] %>%
#    map(~plot_gene(.x, dds, x_var="mutation", color_by="treatment", show_plot = FALSE)$plot)
#  save_plots_to_pdf(plots, file="genes.pdf", ncol=5, nrow=5)

## ---- eval=FALSE--------------------------------------------------------------
#  # design variables need to be factors
#  # since we update the design of the dds
#  # object, we have to do it manually
#  dds$mutation <- as.factor(dds$mutation)
#  dds$treatment <- as.factor(dds$treatment)
#  design(dds) <- ~ mutation + treatment
#  
#  dds <- DESeq(dds, parallel = T)
#  plotDispEsts(dds)

## ---- eval=FALSE--------------------------------------------------------------
#  # get testing results
#  de_res <- lfcShrink(dds, coef = "mutation_WT_vs_D538G", lfcThreshold = log2(1.5), type = "normal", parallel = TRUE)
#  
#  # MA plot
#  plot_ma(de_res, dds, annotate_top_n = 5)

## ---- echo=FALSE--------------------------------------------------------------
plot_ma(T47D_diff_testing, dds, annotate_top_n = 5)

## ---- eval=FALSE--------------------------------------------------------------
#  plot_ma(de_res, dds, highlight_genes = c("CLEC2B", "PAGE5", "GAPDH"))

## ---- echo=FALSE--------------------------------------------------------------
plot_ma(T47D_diff_testing, dds, highlight_genes = c("CLEC2B", "PAGE5", "GAPDH"))

## -----------------------------------------------------------------------------
sessionInfo()

