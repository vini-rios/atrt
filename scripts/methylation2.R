library(GEOquery)
library(tidyverse)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(DMRcate)

# Downloading and unziping data
getGEOSuppFiles("GSE70460")
files_dir <- "GSE70460/files/"

untar("GSE70460/GSE70460_RAW.tar", exdir = files_dir)
idat_files <- list.files(files_dir, pattern = "idat.gz$", full = TRUE)
sapply(idat_files, R.utils::gunzip, overwrite = TRUE)

# From the article ###
### "With the minfi package the same preprocessing steps as in
### Illumina's Genomestudio software were performed."
raw <- read.metharray.exp(base = "GSE70460/files")

mset <- preprocessIllumina(raw, bg.correct = TRUE, normalize = "controls")

### "Additionally, a correction for batch effects was performed.
### Batch effects were estimated by fitting a linear model to the
### log2 transformed intensity values of the methylated and unmethylated channel.
### After removing the component due to the batch effect the residuals were
### back transformed to intensity scale and methylation beta values were calculated
### as described in Illumina’s protocols"

# Here we utlized limma for batch effect removal using the Sentrix ID.

meth <- getMeth(mset)
unmeth <- getUnmeth(mset)


log_meth <- log2(meth + 1)
log_unmeth <- log2(unmeth + 1)


idat_files <- list.files(files_dir, pattern = "_Grn.idat$")

batch <- sub(".*_([0-9]+)_R[0-9]+C[0-9]+.*", "\\1", idat_files)
chip_position <- sub(".*_([0-9]+_R[0-9]+C[0-9]+).*", "\\1", idat_files)

# 63 uniques Sentrix Ids
length(unique(batch))

log_meth_corrected <- removeBatchEffect(log_meth, batch = batch)
log_unmeth_corrected <- removeBatchEffect(log_unmeth, batch = batch)

meth_corrected <- 2^log_meth_corrected - 1
unmeth_corrected <- 2^log_unmeth_corrected - 1

beta <- meth_corrected / (meth_corrected + unmeth_corrected + 100)


### "All samples were checked for expected and unexpected genotype matches
### by pairwise correlation of the 65 genotyping probes on the 450k array."
snp_probes <- getSnpBeta(raw)

cor_matrix <- cor(snp_probes, use = "pairwise.complete.obs")

### Removal of probes targeting the X and Y chromosomes (n = 11,551)
anno <- getAnnotation(mset)

keep_probes <- !(anno$chr %in% c("chrX", "chrY"))

sum(keep_probes) - length(keep_probes)

# 11648 probes were removed from our analysis, 97 more than the original.
# The difference might be due to different annotation versions.
packageVersion("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# ‘0.6.1’

rownames(beta)


methylation_clean <- rmSNPandCH(beta,
  dist = 5,
  mafcut = 0.01,
  and = TRUE,
  rmcrosshyb = TRUE,
  rmXY = TRUE
)

probe_vars <- apply(methylation_clean, 1, var, na.rm = TRUE)
top_probes <- names(sort(probe_vars, decreasing = TRUE)[1:7500])
beta_top <- methylation_clean[top_probes, ]

# Hierarchical clustering
# Euclidean distance for probes, 1-Pearson for samples
sample_dist <- as.dist(1 - cor(beta_top, method = "pearson"))
hc <- hclust(sample_dist, method = "average")

# Plot dendrogram
plot(as.dendrogram(hc, hang = 1, label = NULL), ylim = c(0.2,1), yaxt = "n"  )


gse <- getGEO("GSE70460", GSEMatrix = TRUE)
pheno_geo <- pData(gse[[1]])


library(minfi)
library(limma)
library(GEOquery)

# 1. Get metadata and extract Sentrix info
gse <- getGEO("GSE70460", GSEMatrix = TRUE)
pheno_geo <- pData(gse[[1]])

# Extract Sentrix IDs from supplementary files
extract_sentrix <- function(filename) {
  base <- basename(filename)
  parts <- strsplit(base, "_")[[1]]
  return(parts[2])  # Sentrix ID
}

pheno_geo$sentrix_id <- sapply(pheno_geo$supplementary_file, extract_sentrix)

# 2. Match to your samples
sample_names <- sampleNames(mset)
gsm_ids <- sapply(strsplit(sample_names, "_"), `[`, 1)
pheno_matched <- pheno_geo[match(gsm_ids, pheno_geo$geo_accession), ]

# 3. Define batch
batch <- pheno_matched$sentrix_id

# 4. Check batch distribution
cat("Batch summary:\n")
print(table(batch))

# 5. Batch correction
M <- getMeth(mset)
U <- getUnmeth(mset)

log_M <- log2(M + 1)
log_U <- log2(U + 1)

log_M_corrected <- removeBatchEffect(log_M, batch = batch)
log_U_corrected <- removeBatchEffect(log_U, batch = batch)

M_corrected <- 2^log_M_corrected - 1
U_corrected <- 2^log_U_corrected - 1

# Handle negative values from back-transformation
M_corrected[M_corrected < 0] <- 0
U_corrected[U_corrected < 0] <- 0

# 6. Calculate beta values
beta <- M_corrected / (M_corrected + U_corrected + 100)