library(corpcor)
library(openxlsx)
library(tensorflow)
library(keras)
library(data.table)

setwd("/dir")

rm(list=ls()) 

### Loading trait matrix
Y <- read.xlsx("../Supplementary Table 3_Trait_Matrix.xlsx", rowNames = TRUE)

### CV A Indica-Indica
Y <- Y[Y$Population == "Indica",  -c(1)]
target_accessions <- rownames(Y)
for(i in 1:length(target_accessions)){
  target_accessions[i] <- paste(cbind(target_accessions[i], target_accessions[i]), collapse = "_")
}

### Defining the number of epoch and units
units_M <- 50
epochs_M <- 50

# Reading the markers
GRM <- fread("../kinship.txt")
GRM <- as.matrix(GRM)
lines <- GRM[, c(1)]
row.names(GRM) <- lines
GRM <- GRM[, -c(1)]
colnames(GRM) <- lines

GRM <- GRM[target_accessions, target_accessions]
LG <- t(chol(GRM))
Z1 <- model.matrix(~0+as.factor(target_accessions))
Z1G <- Z1%*%LG
X <- Z1G

### Cross validation scheme
### Define folds
k_folds = 5

# Create folds of equal size
set.seed(1)
folds <- cut(seq(1, nrow(Y)), breaks=k_folds, labels=FALSE)

# Retaining true test values to calculate prediction accuracy
true_fold_test_list <- list()

# List to allocate predictions for the test sets of each fold
prediction_list <- list()

for(i in 1:k_folds){
  ### Getting the test and data partitions for this particular fold
  testIndexes <- which(folds==i, arr.ind=TRUE)
  Y.TST <- Y[testIndexes, ]
  Y.TRN <- Y[-testIndexes, ]
  X.TST <- X[testIndexes,]
  X.TRN <- X[-testIndexes,]
  
  true_fold_test_list[[i]] <- Y.TST
  
  n_accessions <- dim(X)[1]
  
  # Add input layer
  input <- layer_input(shape=dim(X.TRN)[2], name="covars")
  
  # Add hidden layers
  base_model <- input %>%
    
    layer_dense(units =units_M, activation="relu") %>%
    
    layer_dropout(rate = 0.1) %>%
    
    layer_dense(units = units_M, activation = "relu") %>%
    
    layer_dropout(rate = 0.1) %>%
    
    layer_dense(units = units_M, activation = "relu") %>%
    
    layer_dropout(rate = 0.1) %>%
    
    layer_dense(units = units_M, activation = "relu") %>%
    
    layer_dropout(rate = 0.1) %>%
    
    layer_dense(units = units_M, activation = "relu") %>%
    
    layer_dropout(rate = 0.1) %>%
    
    layer_dense(units = units_M, activation = "relu") %>%
    
    layer_dropout(rate = 0.1) %>%
    
    layer_dense(units = units_M, activation = "relu") %>%
    
    layer_dropout(rate = 0.1) %>%
    
    layer_dense(units = units_M, activation = "relu") %>%
    
    layer_dropout(rate = 0.1) %>%
    
    layer_dense(units = units_M, activation = "relu") %>%
    
    layer_dropout(rate = 0.1) %>%
    
    layer_dense(units = units_M, activation = "relu") %>%
    
    layer_dropout(rate = 0.1)
  
  # Add output for trait 1
  
  yhat1 <- base_model %>%
    
    layer_dense(units = 1, name="yhat1")
  
  # Add output 2
  
  yhat2 <- base_model %>%
    
    layer_dense(units = 1, name="yhat2")
  
  # Add output 3
  
  yhat3 <- base_model %>%
    
    layer_dense(units = 1, name="yhat3")
  
  # Add output 4
  
  yhat4 <- base_model %>%
    
    layer_dense(units = 1, name="yhat4")
  
  # Add output 5
  
  yhat5 <- base_model %>%
    
    layer_dense(units = 1, name="yhat5")
  
  # Add output 6
  
  yhat6 <- base_model %>%
    
    layer_dense(units = 1, name="yhat6")
  
  # Add output 7
  
  yhat7 <- base_model %>%
    
    layer_dense(units = 1, name="yhat7")
  
  # Add output 8
  
  yhat8 <- base_model %>%
    
    layer_dense(units = 1, name="yhat8")
  
  # Add output 9
  
  yhat9 <- base_model %>%
    
    layer_dense(units = 1, name="yhat9")
  
  # Add output 10
  
  yhat10 <- base_model %>%
    
    layer_dense(units = 1, name="yhat10")
  
  # Build multi-output model
  
  model <- keras_model(input,list(yhat1,yhat2,yhat3, yhat4, yhat5, yhat6, 
                                  yhat7, yhat8, yhat9, yhat10)) %>%
    
    compile(optimizer = "rmsprop",
            
            loss="mse",
            
            metrics="mae",
            
            loss_weights=c(0.1,0.1,0.1, 0.1, 0.1, 0.1, 0.1,
                           0.1, 0.1 ,0.1))
  
  # fit model
  model_fit <- model %>%
    
    fit(x=X.TRN,
        y=list(Y.TRN[,1],Y.TRN[,2],Y.TRN[,3],Y.TRN[,4],Y.TRN[,5],Y.TRN[,6],
               Y.TRN[,7],Y.TRN[,8],Y.TRN[,9], Y.TRN[,10]),
        epochs=epochs_M,
        batch_size = 50,
        verbose=TRUE)
  
  # Predicting values for the test set
  Yhat <- predict(model, X.TST)
  
  prediction_list[[i]] <- Yhat
}

saveRDS(prediction_list, file="prediction_deeplearning_A_Indica_Indica.RData")

saveRDS(true_fold_test_list, file="true_fold_test_list_deeplearning_A_Indica_Indica.RData")
