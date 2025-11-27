library(GEOquery)
library(tidyverse)
library(readxl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(DMRcate)


methylation_data <- getGEO("GSE70460")
methylation_meta <- getGEO("GSE70460", GSEMatrix = FALSE)

# dkfz to gsm using GEO
dkfz_to_gsm <- purrr::map(methylation_meta@gsms,
    function(x) {
      x@header$title
    }
  ) %>%
  stack() %>%
  dplyr::rename(sample_name = values) %>%
  mutate(dkfz_number = str_extract(sample_name, "(?<=_)[^_]+$"),
    dkfz_number = as.numeric(dkfz_number))

methylation_data <- methylation_data[[1]]

methylation_df <- exprs(methylation_data)




grset <- makeGenomicRatioSetFromMatrix(
  mat = methylation_df,
  array = "IlluminaHumanMethylation450k",
  annotation = "ilmn12.hg19",
  what = "Beta"
)

annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

methylation_clean <- rmSNPandCH(methylation_df,
  dist = 5,
  mafcut = 0.01,
  and = TRUE,
  rmcrosshyb = TRUE,
  rmXY = TRUE
)















annotation_probes <- annotation[rownames(methylation_df),]

keep <- rep(TRUE, nrow(methylation_df))
sum(!keep)
keep <- keep & !(annotation_probes$chr %in% c("chrX", "chrY"))
# Removal of probes targeting the X and Y chromosomes (n = 11,551)
sum(!keep)
# We removed 97 more probes

beta_filtered <- dropLociWithSnps(grset, snps = c("CpG", "SBE"), 
                                   maf = 0.01)

granges(beta_filtered)
granges(grset)




sum(!keep)

methylation_clean <- methylation_df[complete.cases(methylation_df), ]

# Unsupervised hierarchical clustering was performed using the 10,000
# most variable CG-probes over the whole cohort.

probe_vars <- apply(methylation_clean, 1, var)

top_probes <- names(sort(probe_vars, decreasing = TRUE)[1:7500])

methylation_top <- methylation_clean[top_probes, ]

# Clustering steps
# Distance measure was ... 1-Pearson for the samples, average linkage.
sample_dist <- as.dist(1 - cor(methylation_top, method = "pearson"))

hc <- hclust(sample_dist, method = "average")
plot(as.dendrogram(hc, hang = 1, label = NULL), ylim = c(0.2,1), yaxt = "n"  )
rect.hclust(hc, k = 3, border = "red")
my_clusters <- cutree(hc, k = 3)

# Following steps are required to compare out results to the original 

# Creating comparison table matching dkfz from geo to mmc2.xlsx
sup_data <- read_excel("data/mmc2.xlsx", sheet = 2)[2:194, ]
colnames(sup_data) <- sup_data[1,]
sup_data <- sup_data[2:193,] %>%
  janitor::clean_names() %>%
  left_join(., dkfz_to_gsm, by = "sample_name")

sample_order <- colnames(methylation_top)




comparison <- data.frame(
  sample = sample_order,
  my_cluster = my_clusters,
  stringsAsFactors = FALSE
)

comparison$original_subgroup <- sup_data$subgrouping_based_on_450k_methylation_data[
  match(comparison$sample, sup_data$ind)
]

# Checking how many samples matched
sum(!is.na(comparison$original_subgroup))

comparison$my_cluster_group <- case_when(
  comparison$my_cluster == 1 ~ "SHH",
  comparison$my_cluster == 2 ~ "MYC",
  comparison$my_cluster == 3 ~ "TYR"

)

comparison$comparison <- comparison$original_subgroup == comparison$my_cluster_group

sum(comparison$comparison)

comparison[comparison$comparison == FALSE,]

comp <- comparison %>%
  rename(ind = sample) %>%
  left_join(., dkfz_to_gsm, by = "ind")

comparison %>%
  group_by(original_subgroup) %>%
  summarise(count = count(sample))
 