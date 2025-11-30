rule get_microarray_dataset:
    output: "data/processed/microarray/mas5/pascal.csv",
            "data/processed/microarray/mas5/birks.csv",
            "data/processed/microarray/mas5/wang.csv",
            "data/processed/microarray/mas5/amani.csv",
    conda: "../../workflow/envs/microarray.yml"
    script: "../../scripts/get_microarray_dataset.R"

rule get_microarray_metadata:
    output: "data/raw/microarray/metadata/wang_metadata.csv",
            "data/raw/microarray/metadata/amani_metadata.csv",
            "data/raw/microarray/metadata/original_geo_metadata.csv",
            "data/raw/microarray/metadata/mmc2_metadata.csv",

    conda: "../../workflow/envs/microarray.yml"
    script: "../../scripts/get_microarray_metadata.R"

rule batch_correct:
    input:  "data/processed/microarray/mas5/pascal.csv",
            "data/processed/microarray/mas5/birks.csv",
            "data/processed/microarray/mas5/wang.csv",
            "data/processed/microarray/mas5/amani.csv",

    output: "data/processed/microarray/batch/original_study/var_collapsed.rds",
            "data/processed/microarray/batch/original_study/mean_collapsed.rds",
            "data/processed/microarray/batch/original_study/max_collapsed.rds",
            "data/processed/microarray/batch/original_study/test_var_collapsed.rds",
            "data/processed/microarray/batch/original_study/test_mean_collapsed.rds",
            "data/processed/microarray/batch/original_study/test_max_collapsed.rds",

    conda: "../../workflow/envs/microarray.yml"

    script: "../../scripts/batch_correct_2_steps_original.R" 

rule reproduce_original_study:
    input:  "data/raw/microarray/metadata/original_geo_metadata.csv",
            "data/processed/microarray/mas5/pascal.csv",
            "data/processed/microarray/mas5/birks.csv",
            "data/processed/microarray/batch/original_study/var_collapsed.rds",
            "data/processed/microarray/batch/original_study/mean_collapsed.rds",
            "data/processed/microarray/batch/original_study/max_collapsed.rds",

    output: "results/top_1500_var_genes.csv",
            "results/metadata_clusters.csv",
            "results/degs_original_f20.csv",

    conda: "../../workflow/envs/microarray.yml"

    script: "../../scripts/original_study_analysis.R"

rule train_random_forest:
    input:  "results/degs_original_f20.csv",
            "results/metadata_clusters.csv",
            "data/processed/microarray/batch/original_study/var_collapsed.rds",
            "results/top_1500_var_genes.csv",

    output: "models/trained/random_forest.rds",
            "models/trained/top1500_var/random_forest.rds"

    conda: "../../workflow/envs/microarray.yml"

    script: "../../scripts/train_model/train_random_forest.R"

rule train_xgboost:
    input:  "results/degs_original_f20.csv",
            "results/metadata_clusters.csv",
            "data/processed/microarray/batch/original_study/var_collapsed.rds",
            "results/top_1500_var_genes.csv",

    output: "models/trained/xgboost.rds",
            "models/trained/top1500_var/xgboost.rds",

    conda: "../../workflow/envs/microarray.yml"

    script: "../../scripts/train_model/train_xgboost.R"

rule train_plsda:
    input:  "results/degs_original_f20.csv",
            "results/metadata_clusters.csv",
            "data/processed/microarray/batch/original_study/var_collapsed.rds",
            "results/top_1500_var_genes.csv",

    output: "models/trained/plsda.rds",
            "models/trained/top1500_var/plsda.rds",

    conda: "../../workflow/envs/microarray.yml"

    script: "../../scripts/train_model/train_plsda.R"

rule train_glmnet:
    input:  "results/degs_original_f20.csv",
            "results/metadata_clusters.csv",
            "data/processed/microarray/batch/original_study/var_collapsed.rds",
            "results/top_1500_var_genes.csv",

    output: "models/trained/glmnet.rds",
            "models/trained/top1500_var/glmnet.rds",
    
    conda: "../../workflow/envs/microarray.yml"

    script: "../../scripts/train_model/train_glmnet.R"

rule predict_degs:
    input:  "models/trained/random_forest.rds",
            "models/trained/xgboost.rds",
            "models/trained/plsda.rds",
            "models/trained/glmnet.rds",
            "results/metadata_clusters.csv",
            "results/degs_original_f20.csv",
            "data/processed/microarray/batch/original_study/test_var_collapsed.rds",

    output: "results/prediction_degs.csv"

    conda: "../../workflow/envs/microarray.yml"

    script: "../../scripts/predict/predict_degs.R"

rule predict_top_var:
    input:  "models/trained/top1500_var/random_forest.rds",
            "models/trained/top1500_var/xgboost.rds",
            "models/trained/top1500_var/plsda.rds",
            "models/trained/top1500_var/glmnet.rds",
            "results/metadata_clusters.csv",
            "results/top_1500_var_genes.csv",
            "data/processed/microarray/batch/original_study/test_var_collapsed.rds",

    output: "results/prediction_top_var.csv"
    
    conda: "../../workflow/envs/microarray.yml"

    script: "../../scripts/predict/predict_top_var.R"