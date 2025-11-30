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
        expand(
            os.path.join(OUTPUT_DIR, "raw", "{project}_download_complete.flag"),
            project=PROJECT,
        ),
        expand(
            os.path.join(OUTPUT_DIR, "raw", "{project}_all_files_verified.flag"),
            project=PROJECT,
        ),
        expand("results/fastp/{project}_fastp_complete.flag", project=PROJECT),
        "references/homo_sapiens.gtf.gz",
        [
            target
            for project in PROJECT
            for target in get_all_quantification_outputs(project)
        ],
        expand(os.path.join(OUTPUT_DIR, "{project}.metadata.clean"), project=PROJECT),


