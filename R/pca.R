#' Plot results of a principal component analysis
#'
#' @param obj A (features x samples) matrix or SummarizedExperiment object
#' @param PC_x The PC to show on the x-axis.
#' @param PC_y The PC to show on the y-axis.
#' @param n_feats Number of top-variable features to include.
#' @param scale_feats Whether to scale the features.
#' @param na_frac Only consider features with the stated maximum fraction of NAs. NAs will be mean-imputed for PCA.
#' @param metadata A data.frame used for annotating samples. `rownames(metadata)` must match `colnames(obj)`.
#' @param color_by Variable by which to color points. Must be a column in metadata or in `colData(obj)`.
#' @param shape_by Variable by which to color points. Must be a column in metadata or in `colData(obj)`.
#' @param point_alpha alpha value of `geom_point()`
#' @param point_rel_size relative size of `geom_point()`
#' @param show_plot Whether to show the plot or not
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(100*6), ncol=6)
#' data <- t(t(data)+c(-1, -1.1, -1.2, 1, 1.1, 1.2))
#' plot_pca(data)
#'
#' @return The function displays the plot and returns invisible a list of the plot,
#' the data.frame to make the plot, the vector of percentages of variance explained and the loadings matrix.
#' @export
plot_pca <- function(obj, PC_x = 1, PC_y = 2, n_feats = 500, scale_feats = FALSE, na_frac = .3, metadata = NULL, color_by = NULL, shape_by = NULL, point_alpha = .7, point_rel_size = 2, show_plot = TRUE) {

  # prevent 'no visible binding for global variable' package warnings
  sample_id <- feature_id <- NULL

  if (inherits(obj, "SummarizedExperiment")) {
    mat <- tryCatch(as.matrix(assay(obj)), error = function(e) {
      stop("Assay cannot be converted to a matrix.")
    })
    col_data <- as_tibble(colData(obj))
    row_data <- as_tibble(rowData(obj), rownames = "feature_id")
  } else if (inherits(obj, "matrix")) {
    mat <- obj
    col_data <- NULL
    row_data <- NULL
  } else {
    stop("Object is not a matrix or SummarizedExperiment")
  }

  # subset to features with maximum specified fraction of NAs
  mat <- mat[rowSums(is.na(mat)) / ncol(mat) <= na_frac, ]
  # subset to top-variable features
  mat <- mat[matrixStats::rowVars(mat, na.rm = T) %>%
    order(decreasing = T) %>%
    head(n_feats), ]
  # center rows
  mat_centered <- mat %>%
    t() %>%
    scale(scale = scale_feats) %>%
    t()
  # mean impute NAs
  mat_centered[is.na(mat_centered)] <- 0
  # do PCA
  pca <- prcomp(t(mat_centered))
  pca_data <- as_tibble(pca$x, rownames = "sample_id")
  var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  # join annotation
  if (!is.null(col_data)) {
    pca_data <- bind_cols(pca_data, suppressWarnings(select(col_data, -tidyselect::any_of("sample_id"))))
  }
  if (!is.null(metadata)) {
    pca_data <- bind_cols(pca_data, suppressWarnings(select(metadata, -tidyselect::any_of("sample_id"))))
  }

  # first annotation, then PCs
  pca_data <- pca_data %>%
    select(sample_id, !tidyselect::starts_with("PC"), tidyselect::starts_with("PC"))

  loadings_data <- as_tibble(pca$rotation, rownames = "feature_id")
  # add row_data columns from SE object
  if (!is.null(row_data)) {
    loadings_data <- left_join(loadings_data, row_data, by = "feature_id") %>%
      select(feature_id, !tidyselect::starts_with("PC"), tidyselect::starts_with("PC"))
  }

  plot <- pca_data %>%
    ggplot(aes_string(paste0("PC", PC_x), paste0("PC", PC_y), label = "sample_id", color = color_by, shape = shape_by)) +
    geom_point(size = rel(point_rel_size), alpha = point_alpha) +
    coord_fixed() +
    labs(x = paste0("PC", PC_x, " (", round(var_exp[PC_x],1), "%)"), y = paste0("PC", PC_y, " (", round(var_exp[PC_y],1), "%)")) +
    cowplot::theme_cowplot()

  if (show_plot) {
    print(plot)
  }

  invisible(list("plot" = plot, "data" = pca_data, "var_exp" = var_exp, "loadings" = loadings_data))
}

#' Plot loadings of a principal component
#' @param pca_res A result returned from `plot_pca()`
#' @param PC Number of the principal component to plot
#' @param color_by Variable (column in `pca_res$loadings`) to color points by.
#' @param annotate_top_n Annotate the top n features with positive or negative loading
#' @param highlight_genes Vector of gene names or gene IDs to highlight on the plot (overwrites top_n annotation)
#' @param show_plot Whether to show the plot
#'
#' @return The function displays the loadings plot and returns invisible a list of the plot,
#' the data.frame of the PCA loadings.
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(100*6), ncol=6)
#' data <- t(t(data)+c(-1, -1.1, -1.2, 1, 1.1, 1.2))
#' pca_res <- plot_pca(data)
#' plot_loadings(pca_res)
#' @export
plot_loadings <- function(pca_res, PC = 1, color_by = NULL, annotate_top_n = 0, highlight_genes = NULL, show_plot = TRUE) {

  # prevent 'no visible binding for global variable' package warnings
  feature_id <- ylim <- gene_name <- NULL

  loadings_df <- pca_res$loadings %>%
    select(!tidyselect::starts_with("PC"), PC = paste0("PC", PC)) %>%
    mutate(PC_rank = rank(PC))

  # top n features by loading
  annotate_top_n_features <- loadings_df %>%
    filter(
      (rank(PC) <= annotate_top_n) |
        (rank(-PC) <= annotate_top_n)
    ) %>%
    pull(feature_id)

  if (is.null(highlight_genes)) {
    plot <- loadings_df %>%
      ggplot(aes_string(x = "PC_rank", y = "PC", label = "gene_name", color = color_by)) +
      geom_point() +
      # highlight top n
      ggrepel::geom_text_repel(data = loadings_df %>% filter(
        feature_id %in% annotate_top_n_features
      )) +
      labs(x = paste0("rank(PC", PC, " loading)"), y = paste0("PC", PC, " loading")) +
      ylim(c(-max(abs(loadings_df$PC)), max(abs(loadings_df$PC)))) +
      cowplot::theme_cowplot()
  } else {
    plot <- loadings_df %>%
      ggplot(aes_string(x = "PC_rank", y = "PC", label = "gene_name", color = color_by)) +
      geom_point() +
      gghighlight::gghighlight(gene_name %in% highlight_genes | feature_id %in% highlight_genes,
        label_key = gene_name, use_group_by = F,
        label_params = list(fill = NA, label.size = NA, fontface = "bold")
      ) +
      labs(x = paste0("rank(PC", PC, " loading)"), y = paste0("PC", PC, " loading")) +
      ylim(c(-max(abs(loadings_df$PC)), max(abs(loadings_df$PC)))) +
      cowplot::theme_cowplot()
  }

  if (show_plot) {
    print(plot)
  }

  invisible(list("plot" = plot, "data" = loadings_df))
}
