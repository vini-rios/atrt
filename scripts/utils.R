library(tidyverse)
library(GEOquery)
library(matrixStats)

annot_geo <- function(data, gse = NULL, gpl = NULL) {
  # Example gse = "GSE28026"
  # Example gpl = "GPL570"
  if (!is.null(gse) && !is.null(gpl)) {
    stop("Provide gse or gpl, but not both")
  } else if (!is.null(gse)) {
    gset <- getGEO(gse, GSEMatrix = FALSE, getGPL = TRUE)
    gpl_id <- names(gset@gpls)
    gpl_anot <- getGEO(gpl_id, AnnotGPL = TRUE)

    annot_table <- gpl_anot@dataTable@table

    output <- annot_table %>%
      mutate(
        Gene = sapply(`Gene symbol`, function(x) {
          strsplit(x, split = "///")[[1]][1]
        }),
        Gene_alt = sapply(`Gene symbol`, function(x) {
          strsplit(x, split = "///")[[1]][2]
        })
      ) %>%
      select("ID", "Gene", "Gene_alt") %>%
      right_join(., data, by = "ID")

    output
  } else if (!is.null(gpl)) {
    gpl_anot <- getGEO(gpl, AnnotGPL = TRUE)
    annot_table <- gpl_anot@dataTable@table

    output <- annot_table %>%
      mutate(
        Gene = sapply(`Gene symbol`, function(x) {
          strsplit(x, split = "///")[[1]][1]
        }),
        Gene_alt = sapply(`Gene symbol`, function(x) {
          strsplit(x, split = "///")[[1]][2]
        })
      ) %>%
      select("ID", "Gene", "Gene_alt") %>%
      right_join(., data, by = "ID")

    output
  } else {
    stop("Provide gse or gpl")
  }
}

log2_and_matrix <- function(m_data) {
  m_data[, 4:length(m_data)] <- log2(m_data[, 4:length(m_data)] + 1)

  output <- m_data %>%
    column_to_rownames("ID") %>%
    dplyr::select(-Gene, -Gene_alt) %>%
    as.matrix()

  output
}

remove_empty_probes <- function(target, gene_df) {
  gene_corrected <- cbind(gene_df["Gene"], target)
  output <- gene_corrected %>%
    mutate(variance = rowVars(as.matrix(across(where(is.numeric))))) %>%
    arrange(desc(variance)) %>%
    drop_na(Gene) %>%
    select(-variance, -Gene)

  output
}

colapse_probes_gene<- function(target, gse = NULL, gpl = NULL, method = "variance") {
  target_df <- as.data.frame(target) %>%
    tibble::rownames_to_column("ID")
  annot_target <- annot_geo(target_df, gse = gse, gpl = gpl)
  if (method == "variance") {
    collapsed_target <- annot_target %>%
      mutate(variance = rowVars(as.matrix(across(where(is.numeric))))) %>%
      group_by(Gene) %>%
      slice_max(order_by = variance, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(variance)) %>%
      column_to_rownames("Gene") %>%
      dplyr::select(-variance, -ID, -"Gene_alt")
  } else if (method == "mean") {
    collapsed_target <- annot_target %>%
      group_by(Gene) %>%
      summarize(across(where(is.numeric), mean)) %>%
      ungroup() %>%
      column_to_rownames("Gene")
  } else if (method == "max") {
    collapsed_target <- annot_target %>%
      mutate(mean = rowMeans(as.matrix(across(where(is.numeric))))) %>%
      group_by(Gene) %>%
      slice_max(order_by = mean, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      column_to_rownames("Gene") %>%
      dplyr::select(-mean, -ID, -"Gene_alt")
  } else {
    stop("method not recognized")
  }


  collapsed_target
}