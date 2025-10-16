#!/usr/bin/env Rscript
library(tidyverse)
# If ran manually do so from the project root directory (atrt directory) 

arguments <- commandArgs(trailingOnly = TRUE)
project <- arguments[1]

meta_data <- read_tsv(paste0("data/", project, ".metadata")) %>%
  separate(
    description,
    into = c("instrument",
             "name",
             "data"),
    sep = ":"
  ) %>%
  dplyr::select(-instrument, -name)

if (project == "PRJNA565377") {
  meta_data <- meta_data %>%
    separate(
      data,
      into = c("drop1",
               "treatment",
               "drop2",
               "drop3"),
      sep = ","
    ) %>%
    dplyr::select(-drop1, -drop2, -drop3) %>%
    mutate(project = "PRJNA565377")

  write_tsv(meta_data, file = "data/PRJNA565377.metadata.clean")

} else if (project == "PRJNA427332") {
  meta_data <- meta_data %>%
    separate(
      data,
      into = c("drop1", "data"),
      sep = " "
    ) %>%
    separate(
      data,
      into = c("treatment", "drop2", "sample"),
      sep = "_"
    ) %>%
    dplyr::select(-drop1, -drop2) %>%
    mutate(project = "PRJNA427332") %>%
    mutate(treatment = recode(treatment,
                              "nkx" = "NKX2.5",
                              "nkx-31" = "NKX2.5_and_CD31",
                              "31" = "CD31")
    )

  write_tsv(meta_data, file = "data/PRJNA427332.metadata.clean")
}
