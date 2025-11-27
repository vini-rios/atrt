library(tximport)
library(org.Hs.eg.db)


files_path <- list.files("results/quantification/salmon/PRJNA565377",
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
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = TRUE
)
count_matrix <- txi$counts

rownames(count_matrix) <- mapIds(org.Hs.eg.db,
      keys = rownames(count_matrix),
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

write.csv(count_matrix, file = "data/rnaseq/rnaseq.csv")
