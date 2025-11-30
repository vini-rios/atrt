library(randomForest)
library(xgboost)
library(pls)
library(glmnet)
library(tidyverse)
library(ggplot2)

# Load models
random_forest <- readRDS("models/trained/top1500_var/random_forest.rds")
xgboost <- readRDS("models/trained/top1500_var/xgboost.rds")
plsda <- readRDS("models/trained/top1500_var/plsda.rds")
glmnet <- readRDS("models/trained/top1500_var/glmnet.rds")

# Loading cluster information from original_study_analysis.R
metadata <- read.csv("results/metadata_clusters.csv", row.names = 1)
metadata$clusters <- factor(metadata$clusters)

# Loading Top 1500 most variance (responsible for finding the 3 subgroups)
top_var <- read.csv("results/top_1500_var_genes.csv")$x

# Load test dataset
test_data <- readRDS("data/processed/microarray/batch/original_study/test_var_collapsed.rds")

test_df <- test_data[top_var,]


# Random Forest
predictions_rf <- predict(random_forest, t(test_df))

pred_prob_rf <- predict(random_forest, t(test_df), type = "prob")

confidence_rf <- apply(pred_prob_rf, 1, max)

results_rf <- data.frame(
  sample = colnames(test_df),
  technology = "microarray",
  random_forest_cluster = predictions_rf,
  random_forest_confidence = confidence_rf
)
rownames(results_rf) <- NULL

# XGBoots
cluster_names <- levels(metadata$clusters)
xgb_test <- xgb.DMatrix(data = t(test_df))

pred_prob_xgb <- predict(xgboost, xgb_test, type = "prob")

confidence_xgb <- apply(pred_prob_xgb, 1, max)

xgb_results <- data.frame(
  sample = colnames(test_df),
  xgboost_cluster = cluster_names[apply(pred_prob_xgb, 1, which.max)],
  xgboost_confidence = confidence_xgb
)

# PLS-DA (Partial Least-Squares Discriminant Analysis)
pls_test <- t(test_df)

pred_scores_pls <- predict(plsda, newdata = pls_test, ncomp = 3)
test_proj_pls <- predict(plsda, newdata = pls_test, type = "scores")

pred_class_pls <- apply(pred_scores_pls, 1, which.max)

y_train <- metadata$clusters
pred_classes_pls   <- levels(y_train)[pred_class_pls]

get_probs <- function(x) { exp(x) / sum(exp(x)) }
pls_probs <- t(apply(pred_scores_pls, 1, get_probs))

pls_results <- data.frame(
  sample = rownames(pls_test),
  pls_cluster = pred_classes_pls,
  pls_confidence = apply(pls_probs, 1, max)
)

test_scores_pls <- data.frame(
  Comp1 = test_proj_pls[, 1],
  Comp2 = test_proj_pls[, 2],
  Group = pred_classes_pls,
  Set = "Test_data"
)
train_scores <- data.frame(
  Comp1 = scores(plsda)[, 1],
  Comp2 = scores(plsda)[, 2],
  Group = y_train,
  Set   = "Pascal"
)
all_scores <- rbind(test_scores_pls, train_scores)
ggplot(all_scores, aes(x = Comp1, y = Comp2, color = Group, shape = Set)) +
  geom_point(size = 3, alpha = 0.7) +                 
  stat_ellipse(data = train_scores, level = 0.95) +   
  theme_minimal() +
  labs(title = "PLS-DA: Training vs Test Data",
       subtitle = "Ellipses are 95% Confidence Intervals.",
       x = "Latent Component 1", y = "Latent Component 2")

# Penalized Logistic Regression - glmnet
y_train <- metadata$clusters

glmnet_test <- t(test_df)
pred_probs_glmnet <- predict(glmnet, glmnet_test, type = "response", s = "lambda.min")

glmnet_results <- data.frame(
  sample = colnames(test_df),
  glmnet_cluster = colnames(pred_probs_glmnet)[apply(pred_probs_glmnet, 1, which.max)],
  glmnet_confidence = apply(pred_probs_glmnet, 1, max)
)

all(results_rf$sample == xgb_results$sample)
all(results_rf$sample == pls_results$sample)
all(results_rf$sample  == glmnet_results$sample)

overall <- cbind(
  results_rf,
  xgb_results %>% select(-sample),
  pls_results %>% select(-sample),
  glmnet_results %>% select(-sample)
)

write.csv(overall, file = "results/prediction_top_var.csv")