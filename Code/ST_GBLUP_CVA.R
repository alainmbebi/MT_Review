library(BGLR)
library(corpcor)
library(openxlsx)
library(data.table)

rm(list=ls()) 

### Loading trait matrix
Y <- read.xlsx("../Supplementary Table 3_Trait_Matrix.xlsx.xlsx", rowNames = TRUE)

### CV A Indica-Indica
Y <- Y[Y$Population == "Indica",  -c(1)]
target_accessions <- rownames(Y)
for(i in 1:length(target_accessions)){
  target_accessions[i] <- paste(cbind(target_accessions[i], target_accessions[i]), collapse = "_")
}
n_traits <- ncol(Y)

#BGLR parameters
itermax=10000
burning = 5000

# Reading the markers
X <- fread("../Supplementary File 1_markers.txt")
colnames(X)[1] <- "Marker"
X <- X[X$Marker %in% target_accessions, -c(1)]

### Setting up Cross Validation approach 10
### Define folds
k_folds <- 5

# Create folds of equal size
set.seed(1)
folds <- cut(seq(1, nrow(Y)), breaks=k_folds, labels=FALSE)

# Retaining true test values to calculate prediction accuracy
true_fold_test_list <- list()

# List to allocate predictions for the test sets of each fold
prediction_list <- list()

for(i in 1:k_folds){
  
  sprintf("Current fold: %s", i)
  
  ### Getting the test and data partitions for this particular fold
  testIndexes <- which(folds==i,arr.ind=TRUE)
  Y.TST <- Y[testIndexes, ]
  Y.TRN <- Y[-testIndexes, ] 
  X.TST <- X[testIndexes]
  true_fold_test_list[[i]] <- Y.TST
  
  # This will store the results per trait as vector columns
  Y.PRED <- matrix(nrow=nrow(X.TST), ncol=n_traits)
  X.TST <- as.matrix(X.TST)
  
  # ETA
  ETA <- list(list(X=X[-testIndexes,], model='BRR'))
  
  for(l in 1:n_traits){
    sprintf("Current trait for this fold: %s", l)
    fm <- BGLR(y=Y.TRN[,l], ETA=ETA, nIter=itermax, burnIn=burning, verbose=TRUE)
    
    yHat <- fm$mu + as.vector(X.TST%*%fm$ETA[[1]]$b)
    
    Y.PRED[,l] <- yHat
  }
  prediction_list[[i]] <- Y.PRED
}

saveRDS(prediction_list, file="prediction_unitrait_A_Indica_Indica.RData")

saveRDS(true_fold_test_list, file="true_fold_test_list_unitrait_A_Indica_Indica.RData")
