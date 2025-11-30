import os
import glob
import pathlib
import zipfile


configfile: "config.yaml"


include: "workflow/rules/common.smk"
include: "workflow/rules/process_rnaseq_data.smk"
include: "workflow/rules/process_microarray.smk"


PROJECT = ["PRJNA565377"]
OUTPUT_DIR = "data/"


rule all:
    input:
        "results/prediction_top_var.csv",
        "results/prediction_degs.csv",


