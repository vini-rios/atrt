library(tidyverse)
library(affy)
library(GEOquery)


source("scripts/utils.R")

# Functions
# Download raw data from GEO
raw_geo <- function(gse) {
  base_dir <- "data/raw/microarray/"
  files_dir <- paste0(base_dir, gse, "/files")

  getGEOSuppFiles(gse, baseDir = base_dir)

  untar(paste0("data/raw/", gse, "/", gse, "_RAW.tar"), exdir = files_dir)

  raw_data <- ReadAffy(celfile.path = files_dir)

  raw_data
}

# Mas5 normalization
mas5_norm <- function(raw_data) {
  mas5_data <- mas5(raw_data, sc = 100)

  mas5_df <- exprs(mas5_data)

  mas5_df <- rownames_to_column(as.data.frame(mas5_df), var = "ID")

  colnames(mas5_df) <- str_remove(colnames(mas5_df), "_.*")

  mas5_df
}

# Download, mas5 normalization, annotate
process_dataset_mas5 <- function(gse) {
  output <- raw_geo(gse) %>%
    mas5_norm(.) %>%
    annot_geo(., gse = gse)

  output
}

# Here we use GEOquery to download and process microarray data from both
# datasets (GSE70678, GSE28026)
# We also compare the processed raw data to the processed data from GEO
# from GSE70678 to verify our analysis.


# Processed files from GEO 
gset <- getGEO("GSE70678", GSEMatrix = TRUE, getGPL = FALSE)

gset <- gset[[1]]

geo_df  <- exprs(gset)

geo_df <- rownames_to_column(as.data.frame(geo_df), var = "ID")

# Raw files from GEO
raw_data <- raw_geo("GSE70678")

# Normalize data with mas5
mas5_df <- mas5_norm(raw_data)

# Rounding to match processed files from GEO
round_df <- mas5_df
round_df[, 2:50] <- lapply(round_df[, 2:50], function(x) round(x, 1))

# Comparing both = FALSE
all(round_df == geo_df)

# Finding the issue

# First sample and microarray_names have the same mismatches?
matches_sample <- round_df[, 1] == geo_df[, 1]

matches_names <- round_df["ID"] == geo_df["ID"]

miss_match_names <- which(matches_names == FALSE)
miss_match_sample <- which(matches_sample == FALSE)

all(miss_match_names == miss_match_sample)

# Yes. Is this the same pattern across the df?
matches <- round_df == geo_df

all_matches <- rowSums(matches) == dim(matches)[2]

miss_match_all <- which(all_matches == FALSE)

all(miss_match_all == miss_match_names)

length(miss_match_all)


# Yes, 1118 rows are out of order
round_df <- round_df %>%
  arrange(ID)

geo_df <- geo_df %>%
  arrange(ID)

all(round_df == geo_df)

# Sorting both dataframes fixed the issue and both are identical.
# Raw data + mas5 normalization perfectly matches the GEO processed file.

# Annotate table utils.R
pascal <- annot_geo(mas5_df, gse = "GSE70678")

write_csv(pascal, "data/processed/microarray/mas5/pascal.csv")

# Birks "GSE28026",
# Raw data since it was normalized using the gcRMA algorithm instead of mas5
birks_gse <- "GSE28026"

birks <- process_dataset_mas5(birks_gse)

write_csv(birks, "data/processed/microarray/mas5/birks.csv")

# Datasets below were also normalized using different methods
# Wang GSE65132
wang_gse <- "GSE65132"

wang <- process_dataset_mas5(wang_gse)

write_csv(wang, "data/processed/microarray/mas5/wang.csv")

# Amani GSE86574
amani_gse <- "GSE86574"

amani <- process_dataset_mas5(amani_gse)

write_csv(amani, "data/processed/microarray/mas5/amani.csv")
