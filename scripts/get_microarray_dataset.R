library(affy)
library(GEOquery)
library(tidyverse)

# Here we use GEOquery to download and process microarray data from both
# datasets (GSE70678, GSE28026)
# We also compare the processed raw data to the processed data from GEO
# from GSE70678 to verify our analysis.
# Final output is "data/microarray/microarray_dataset.csv"

# Processed files from GEO 
gset <- getGEO("GSE70678", GSEMatrix = TRUE, getGPL = FALSE)

gset <- gset[[1]]

geo_df  <- exprs(gset)

geo_df <- rownames_to_column(as.data.frame(geo_df), var = "ID")

# Raw files from GEO
files_dir <- "GSE70678/files/"

getGEOSuppFiles("GSE70678")

untar("GSE70678/GSE70678_RAW.tar", exdir = files_dir)

raw_data <- ReadAffy(celfile.path = files_dir)

# Normalize data with mas5
mas5_data <- mas5(raw_data, sc = 100)

mas5_df <- exprs(mas5_data)

mas5_df <- rownames_to_column(as.data.frame(mas5_df), var = "ID")

colnames(mas5_df) <- str_remove(colnames(mas5_df), "_.*")

# Rounding to match processed files from GEO
mas5_df[, 2:50] <- lapply(mas5_df[, 2:50], function(x) round(x, 1))

# Comparing both = FALSE
all(mas5_df == geo_df)

# Finding the issue

# First sample and microarray_names have the same mismatches?
matches_sample <- mas5_df[, 1] == geo_df[, 1]

matches_names <- mas5_df["ID"] == geo_df["ID"]

miss_match_names <- which(matches_names == FALSE)
miss_match_sample <- which(matches_sample == FALSE)

all(miss_match_names == miss_match_sample)

# Yes. Is this the same pattern across the df?
matches <- mas5_df == geo_df

all_matches <- rowSums(matches) == dim(matches)[2]

miss_match_all <- which(all_matches == FALSE)

all(miss_match_all == miss_match_names)

length(miss_match_all)


# Yes, 1118 rows are out of order
mas5_df <- mas5_df %>%
  arrange(ID)

geo_df <- geo_df %>%
  arrange(ID)

all(mas5_df == geo_df)

# Sorting both dataframes fixed the issue and both are identical.
# Raw data + mas5 normalization perfectly matches the GEO processed file.

# GPL570 ID to gene symbol
gpl_id <- annotation(gset)
gpl <- getGEO(gpl_id, AnnotGPL = TRUE)

annot_table <- gpl@dataTable@table

# Final df
microarray_df <- annot_table %>%
  mutate(
    Gene = sapply(`Gene symbol`, function(x) {
      strsplit(x, split = "///")[[1]][1]
    }),
    Gene_alt = sapply(`Gene symbol`, function(x) {
      strsplit(x, split = "///")[[1]][2]
    })
  ) %>%
  select("ID", "Gene", "Gene_alt") %>%
  inner_join(., geo_df, by = "ID")



# Birks dataset,
# Raw data since it was normalized using the gcRMA algorithm instead of mas5
files_dir_2 <- "GSE28026/files/"
getGEOSuppFiles("GSE28026")
untar("GSE28026/GSE28026_RAW.tar", exdir = files_dir_2)
raw_birks <- ReadAffy(celfile.path = files_dir_2)

mas5_birks <- mas5(raw_birks, sc = 100)

birks_mas5_df <- exprs(mas5_birks)

birks_mas5_df <- rownames_to_column(as.data.frame(birks_mas5_df), var = "ID")

colnames(birks_mas5_df) <- str_remove(colnames(birks_mas5_df), "_.*")

birks_mas5_df[, 2:19] <- lapply(birks_mas5_df[, 2:19], function(x) round(x, 1))

# Adding birks data set to microarray_df
dim(microarray_df)[1] == dim(birks_mas5_df)[1]
microarray_df <- microarray_df %>%
  inner_join(., birks_mas5_df, by = "ID")

dim(microarray_df)

write_csv(microarray_df, "data/microarray/microarray_dataset.csv")


# Wang GSE65132
files_dir_3 <- "GSE65132/files/"
getGEOSuppFiles("GSE65132")
untar("GSE65132/GSE65132_RAW.tar", exdir = files_dir_3)
raw_wang <- ReadAffy(celfile.path = files_dir_3)
mas5_wang <- mas5(raw_wang, sc = 100)

wang_mas5_df <- exprs(mas5_wang)

wang_mas5_df <- rownames_to_column(as.data.frame(wang_mas5_df), var = "ID")

colnames(wang_mas5_df) <- str_remove(colnames(wang_mas5_df), "_.*")

wang_df <- annot_table %>%
  mutate(
    Gene = sapply(`Gene symbol`, function(x) {
      strsplit(x, split = "///")[[1]][1]
    }),
    Gene_alt = sapply(`Gene symbol`, function(x) {
      strsplit(x, split = "///")[[1]][2]
    })
  ) %>%
  select("ID", "Gene", "Gene_alt") %>%
  inner_join(., wang_mas5_df, by = "ID")


write_csv(wang_df, "data/microarray/wang.csv")

# Amani GSE86574
files_dir_4 <- "GSE86574/files/"
getGEOSuppFiles("GSE86574")
untar("GSE86574/GSE86574_RAW.tar", exdir = files_dir_4)
raw_amani <- ReadAffy(celfile.path = files_dir_4)
mas5_amani <- mas5(raw_amani, sc = 100)

amani_mas5_df <- exprs(mas5_amani)

amani_mas5_df <- rownames_to_column(as.data.frame(amani_mas5_df), var = "ID")

colnames(amani_mas5_df) <- str_remove(colnames(amani_mas5_df), "_.*")

amani_df <- annot_table %>%
  mutate(
    Gene = sapply(`Gene symbol`, function(x) {
      strsplit(x, split = "///")[[1]][1]
    }),
    Gene_alt = sapply(`Gene symbol`, function(x) {
      strsplit(x, split = "///")[[1]][2]
    })
  ) %>%
  select("ID", "Gene", "Gene_alt") %>%
  inner_join(., amani_mas5_df, by = "ID")

write_csv(amani_df, "data/microarray/amani.csv")
