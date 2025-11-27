library(tidyverse)

annot_geo <- function(gse, data) {
  gset <- getGEO(gse, GSEMatrix = FALSE, getGPL = TRUE)
  gpl_id <- names(gset@gpls)
  gpl <- getGEO(gpl_id, AnnotGPL = TRUE)

  annot_table <- gpl@dataTable@table

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
    inner_join(., data, by = "ID")

  output
}

colapse_probes <- function(target, gene_df) {
  gene_corrected <- cbind(gene_df["Gene"], target)
  collapsed_target <- gene_corrected %>%
    mutate(variance = rowVars(as.matrix(across(where(is.numeric))))) %>%
    drop_na(Gene) %>%
    group_by(Gene) %>%
    slice_max(order_by = variance, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(variance)) %>%
    column_to_rownames("Gene") %>%
    dplyr::select(-variance)

  collapsed_target
}
