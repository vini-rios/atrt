library(GEOquery)
library(tidyverse)
library(janitor)
library(readxl)
# Here we create the metadata table
# We use mainly the GEO samples information.

# We need to use mmc2.xlsx
# In order to determine dkfz sample names from Birks dataset (GSE28026)
# And to compare our findings to the original study
# mmc2.xlsx is a suplemental material:
# Table S1. Overview of All ATRTs Used in the Study and Analyses, Related to Figure 1
# From DOI: 10.1016/j.ccell.2016.02.001

# We could have extracted almost all information from mmc2.xlsx
# Only the dkfz_ATRT_# to GSM# couldn't be extracted.
# But we chose to use the information provided by GEO as a learning exercise.


# Creating metadata table
meta_pascal <- getGEO("GSE70678", GSEMatrix = FALSE)
meta_birks <- getGEO("GSE28026", GSEMatrix = FALSE)

get_metadata_geo <- function(gset_meta) {
  metadata <- purrr::map(gset_meta@gsms,
    function(x) {
      x@header$characteristics_ch1
    }
  ) %>%
    stack() %>%
    tidyr::separate(values, into = c("feature", "value"), sep = ": ") %>%
    tidyr::pivot_wider(
      names_from = feature,
      values_from = value
    ) %>%
    janitor::clean_names() 

  dkfz <- purrr::map(gset_meta@gsms,
    function(x) {
      x@header$title
    }
  ) %>%
    stack()

  metadata <- inner_join(metadata, dkfz, by = "ind")

  metadata
}

meta_pascal_tibble <- get_metadata_geo(meta_pascal) %>%
  select("ind", "gender", "tumor_location", "values") %>%
  dplyr::rename(sample_name = values) %>%
  mutate(project = "pascal")

meta_birks_tibble <- get_metadata_geo(meta_birks) %>%
  mutate(project = "birks") %>%
  mutate(tumor_location = case_when(
    tumor_location == "posterior fossa" ~ "infratentorial",
    tumor_location == "frontal lobe" ~ "supratentorial",
    tumor_location == "parietal / occipital lobes" ~ "supratentorial",
    tumor_location == "temporal" ~ "supratentorial",
    TRUE ~ NA
  )
  ) %>%
  mutate(features = paste(
    gender,
    round(as.numeric(age_at_diagnosis_months)/12, digits = 2),
    tumor_location,
    sep = "_"
  )
  ) %>%
  mutate(age_at_diagnosis_years = as.numeric(age_at_diagnosis_months) / 12) %>%
  select(
    ind,
    features,
    gender,
    age_at_diagnosis_months,
    age_at_diagnosis_years,
    tumor_location
  )

# Features is a unique identifier to link dkfz sample name to GSM
length(unique(meta_birks_tibble$features)) == dim(meta_birks_tibble)[1]



# We need the following to determine the dkfz_ATRT_# of the Birks dataset
# and to validate our fiding to the original.

# The dkfz_ATRT_# to GSM# from the 26 samples with
# both methylation data and microarray data are found in the GEO
# So we could determine the 1500 DEGs for clustering without mmc2.xlsx

# Suplemental data from Pascal Table S1 DOI: 10.1016/j.ccell.2016.02.001
# Trying to use download.file() will return HTTP status was '403 Forbidden'
# So the file was manually donwloaded from the link:
# https://www.cell.com/cms/10.1016/j.ccell.2016.02.001/attachment/41d02f67-5661-42b4-aa66-99beeb125458/mmc2.xlsx

sup_data <- read_excel("data/mmc2.xlsx", sheet = 2)[2:194, ]
colnames(sup_data) <- sup_data[1,]
sup_data <- sup_data[2:193,] %>%
  janitor::clean_names()

# Fixed rounding in sup. material to consistently match age in months by birks
# Otherwise we couldn't get a match from GEO and mmc2.xlsx

# Only age in years was changed to the following values
# Had to guess a few values that were not rounding mismatches
# Since age + tumor location + sex was unique this process was doable
# I might be wrong here so I am highlighting the changes
# I also tried to show my reasoning for the changes
# Numbers are out of order because some sample were identified by
# being the last available for tumor location + sex

# *** to highlight biggest changes were I might be wrong
# in parenthesis (age in months/12) and/or comment

# dkfz_ATRT_181 from 1.6 to 1.58 (1.5833)
sup_data$age_at_diagnosis[sup_data$sample_name == "dkfz_ATRT_181"] <- 1.58

# dkfz_ATRT_182 from NA to 0
sup_data$age_at_diagnosis[sup_data$sample_name == "dkfz_ATRT_182"] <- 0

# dkfz_ATRT_186 from 0.6 to 0.58 (0.5833)
sup_data$age_at_diagnosis[sup_data$sample_name == "dkfz_ATRT_186"] <- 0.58

# dkfz_ATRT_190 from 3.12 to 3.17 (3.1667) (mistype?) ***
sup_data$age_at_diagnosis[sup_data$sample_name == "dkfz_ATRT_190"] <- 3.17

# dkfz_ATRT_191 from 1.66 to 1.67 (1.667) (only F supra)
sup_data$age_at_diagnosis[sup_data$sample_name == "dkfz_ATRT_191"] <- 1.67

# dkfz_ATRT_188 from 1.67 to 1.42 (1.4167) (last remaining F) ***
sup_data$age_at_diagnosis[sup_data$sample_name == "dkfz_ATRT_188"] <- 1.42

# dkfz_ATRT_180 from 0.91 to 0.92 (0.9167)
sup_data$age_at_diagnosis[sup_data$sample_name == "dkfz_ATRT_180"] <- 0.92

# dkfz_ATRT_179 from 1.66 to 0.83 (last remaining M infra)
# (Maybe 1.66/2 = 0.83 ?) ***
sup_data$age_at_diagnosis[sup_data$sample_name == "dkfz_ATRT_179"] <- 0.83

# dkfz_ATRT_187 from 1.3 to 2.75 (last remaining M supra)
# (Maybe 1.33*2 = 2.66 ?) ***
sup_data$age_at_diagnosis[sup_data$sample_name == "dkfz_ATRT_187"] <- 2.75

# There are some mismatches in dkfz sample names in mmc2.xlsx and in GEO
sup_material_names <- sup_data %>%
  select(sample_name, subrouping_based_on_affymetrix_gene_expression_data) %>%
  drop_na(subrouping_based_on_affymetrix_gene_expression_data)


# The following sample names are in mmc2.xlsx.
# They have subrouping_based_on_affymetrix_gene_expression_data
# But are not present in GSE70678
sup_material_names$sample_name[!(sup_material_names$sample_name %in% meta_pascal_tibble$sample_name)]

# 25 samples
sum(!(sup_material_names$sample_name %in% meta_pascal_tibble$sample_name))



# Final birks metadata process
meta_birks_final <- sup_data %>%
  filter(previosuly_published_pmid_number == 21946044) %>%
  mutate(age_at_diagnosis = round(as.numeric(age_at_diagnosis), digits = 2)) %>%
  mutate(features = paste(
    gender_f_fe_male_m_male,
    age_at_diagnosis,
    localization_of_primary_tumor,
    sep = "_"
  )
  ) %>%
  left_join(., meta_birks_tibble, by = "features") %>%
  select(
    sample_name,
    ind,
    tumor_location,
    gender,
  ) %>%
  mutate(project = "birks")

# Combining both metadata
metadata <- rbind(meta_pascal_tibble, meta_birks_final) %>%
  dplyr::rename("gsm" = ind) %>%
  left_join(.,
    sup_data %>%
      select(
        sample_name,
        molecular_subgroup_consensus,
        subgrouping_based_on_450k_methylation_data,
        subrouping_based_on_affymetrix_gene_expression_data
      ),
    by = "sample_name"
  )

write_csv(metadata, "data/microarray/microarray_metadata.csv")

# Amani GSE86574
meta_amani <- getGEO("GSE86574", GSEMatrix = FALSE)

meta_amani_tibble <- get_metadata_geo(meta_amani)
