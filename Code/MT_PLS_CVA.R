library(SKM)
library(dplyr)
library(data.table)
library(RhpcBLASctl)
library(openxlsx)
library(Matrix)

rm(list = ls())

# Loading traits
Y <- read.xlsx("../Supplementary Table 3_Trait_Matrix.xlsx", rowNames = TRUE)

### CV A Indica-Indica
Y <- Y[Y$Population == "Indica",  -c(1)]
target_accessions <- rownames(Y)
for(i in 1:length(target_accessions)){
  target_accessions[i] <- paste(cbind(target_accessions[i], target_accessions[i]), collapse = "_")
}

# Reading the markers
GRM <- fread("../kinship.txt")
GRM <- as.matrix(GRM)
lines <- GRM[, c(1)]
row.names(GRM) <- lines
GRM <- GRM[, -c(1)]
colnames(GRM) <- lines

# Transforming markers for the X design matirx
GRM <- GRM[target_accessions, target_accessions]
LG <- t(chol(GRM))
Z1 <- model.matrix(~0+as.factor(target_accessions))
Z1G <- Z1%*%LG
X <- Z1G

############################### Multitrait ####################################
# Folds
k_folds <- 5

colnames(Y)[3] <- "Heading.Date"
colnames(Y)[5] <- "Plant.Height"

# Create folds of equal size
set.seed(1)
folds <- cv_kfold(records_number = nrow(X), k = k_folds)

# Retaining true test values to calculate prediction accuracy
true_fold_test_list <- list()

# List to allocate predictions for the test sets of each fold
prediction_list <- list()

headers <- rownames(X)
colnames(X) <- headers
rownames(Y) <- headers

# Model training and predictions of the ith partition
for(i in seq_along(folds)){
  cat("\n")
  cat("\t\t\t\t*** Fold:", i, " ***\n")
  fold <- folds[[i]]

  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- Y[fold$training,]
  y_testing <- Y[fold$testing,]
  
  true_fold_test_list[[i]] <- y_testing
  
  # Model training
  model <- partial_least_squares(
    x = X_training,
    y = y_training
  )
  
  #Prediction of the testing set
  ncomp_line <- model$optimal_components_num
  predictions <- predict(model, X_testing, components_num = ncomp_line, format = "data.frame")
  
  prediction_list[[i]] <- predictions # For this fold
}
