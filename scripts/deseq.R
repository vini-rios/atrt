library(tidyverse)
library(GenomicFeatures)
library(txdbmaker)
library(tximport)
library(DESeq2)
library(org.Hs.eg.db)
library(RColorBrewer)
library(pheatmap)
library(patchwork)
library(ComplexUpset)

params <- list(project= "PRJNA565377", zero=0)

metadata <- readr::read_tsv(paste0("data/", params$project, ".metadata.clean"))


files_path <- list.files(paste0("results/quantification/salmon/", params$project),
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE)

files_names <- dirname(files_path) %>%
  basename()

names(files_path) <- files_names



gtf <- "references/homo_sapiens.gtf.gz"

txdb <- txdbmaker::makeTxDbFromGFF(gtf)

k <- biomaRt::keys(txdb, keytype = "TXNAME")

tx2gene <- biomaRt::select(txdb,
  keys = k,
  columns = "GENEID",
  keytype = "TXNAME")

colnames(tx2gene) <- c("TXNAME", "GENEID")


txi <- tximport(
  files = files_path,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)


dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = metadata,
  design = ~treatment
)

smallest_group_size <- metadata %>%
    dplyr::count(treatment)%>%
    summarise(min_group_size = min(n)) %>%
    dplyr::pull(min_group_size)


keep <- rowSums(counts(dds) >= 10) >= smallest_group_size
dds <- dds[keep, ]


dds$treatment <- relevel(dds$treatment, ref = "Sol")

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)

FDR_threshold <- 0.05
lfc_threshold <- 1


DESeq2::resultsNames(dds)

comparison = "treatment_MMA_8h_vs_Sol"

results <- DESeq2::results(dds, name = comparison)

results_shrink <- DESeq2::lfcShrink(
    dds,
    coef = comparison,
    type = "apeglm"
)

results_shrink_df <- as.data.frame(results_shrink) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(padj)

results_df <- as.data.frame(results) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(padj)

results_shrink_df$symbol <- mapIds(org.Hs.eg.db,
    keys = results_shrink_df$gene_id,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
)

results_shrink_df$gene_name <- mapIds(org.Hs.eg.db,
    keys = results_shrink_df$gene_id,
    column = "GENENAME",
    keytype = "ENSEMBL",
    multiVals = "first"
)

degs_shrink<- results_shrink_df %>%
    filter(padj < FDR_threshold & abs(log2FoldChange) > lfc_threshold)






plotMA(results, ylim=c(-2,2))
plotMA(results_shrink, ylim=c(-2,2))


volcano_data <- results_df %>%
    mutate(deg = case_when(
      padj <= FDR_threshold & log2FoldChange > lfc_threshold ~ "Upregulated",
      padj <= FDR_threshold & log2FoldChange < -lfc_threshold ~ "Downregulated",
      TRUE ~ "Not-significant"
      )
    ) %>%
    mutate(significant = case_when(
      padj <= FDR_threshold ~"Significant",
      padj > FDR_threshold ~"Non-significant",
      is.na(padj) & !is.na(pvalue) ~ "Low normalized count",
      TRUE ~ "Outlier"
      )
    ) %>%
    mutate(y_display = case_when(
      log2FoldChange > 4 ~ 4,
      log2FoldChange < -4 ~ -4,
      TRUE ~ log2FoldChange
    ),
    shape_type = case_when(
      
      log2FoldChange > 4 ~ "upper",
      log2FoldChange < -4 ~ "lower",
      TRUE ~ "normal"
    )
  )

ma_plot <- ggplot(volcano_data, aes(x = baseMean, y = log2FoldChange, color = significant, alpha = significant)) +
geom_point(data = subset(volcano_data, shape_type == "normal"), 
            size = 2) +
geom_point(data = subset(volcano_data, shape_type == "upper"),
            shape = 17, size = 2.5) +  # upward triangle
geom_point(data = subset(volcano_data, shape_type == "lower"),
            shape = 6, size = 2.5) + 
scale_x_log10() +
theme_minimal() +
labs(x = expression(log[10]~"(Mean of Normalized Counts)"), y = expression("Shrunken"~log[2]~"Fold Change"))+
ylim(c(-5, 5))+
theme(legend.title = element_blank(),
        legend.position = "top")+
scale_color_manual(values = c("Significant" = "red", "Non-significant" = "grey", "Low normalized count" = "#895129", "Outlier" = "#8B008B"))+
scale_alpha_manual(values = c("Significant" = 0.8, "Non-significant" = 0.3, "Low normalized count" = 0.5, "Outlier" = 0.5))+
guides(color = guide_legend(override.aes = list(shape = 16)))



volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = deg, alpha = deg)) +
geom_point(size = 2) +
scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "#6960dbff", 
                                "Not significant" = "gray")) +
scale_alpha_manual(values = c("Upregulated" = 0.8, 
                                "Downregulated" = 0.8, 
                                "Not significant" = 0.5))+
geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed") +
labs(x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"Adjusted P-value")) +
theme_minimal()+
theme(legend.title = element_blank(),
        legend.position = "top")


dds <- DESeq(dds, minReplicatesForReplace=Inf)
res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
results <- res
summary(res)


library(edgeR)
library(limma)

# Get gene-level counts
counts_matrix <- round(txi$counts)

# Create DGEList
dge <- DGEList(counts = counts_matrix, 
               group = metadata$treatment)

# Filter
keep <- filterByExpr(dge, group = metadata$treatment)
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes retained:", sum(keep), "\n")

# TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# Design
design <- model.matrix(~0 + treatment, data = metadata)
colnames(design) <- levels(metadata$treatment)

# Voom
v <- voom(dge, design, plot = TRUE)

# PCA check
plotMDS(v, labels = metadata$treatment)

# Fit
fit <- lmFit(v, design)

# Contrasts
contrasts <- makeContrasts(
  MMA_18h_vs_Sol = MMA_18h - Sol,
  MMA_8h_vs_Sol = MMA_8h - Sol,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Results
res_18h <- topTable(fit2, coef = "MMA_18h_vs_Sol", number = Inf)
res_8h <- topTable(fit2, coef = "MMA_8h_vs_Sol", number = Inf)

cat("\n18h DEGs (FDR < 0.05):", sum(res_18h$adj.P.Val < 0.05), "\n")
cat("18h DEGs (FDR < 0.05, |logFC| > 1):", 
    sum(res_18h$adj.P.Val < 0.05 & abs(res_18h$logFC) > 1), "\n")

cat("\n8h DEGs (FDR < 0.05):", sum(res_8h$adj.P.Val < 0.05), "\n")
cat("8h DEGs (FDR < 0.05, |logFC| > 1):", 
    sum(res_8h$adj.P.Val < 0.05 & abs(res_8h$logFC) > 1), "\n")