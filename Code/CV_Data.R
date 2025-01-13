library(openxlsx)
library(corpcor)
library(BMTME)

rm(list=ls()) 

#### Bringing trait values
preprocess_xlsx <- function(filename, sheet){
  file <- read.xlsx(filename, sheet)
  df <- data.frame(file, stringsAsFactors = False)
  headers <- df[,1]
  df <- df[,-1]
  rownames(df) <- headers
  return(df)
}

traits <- preprocess_xlsx("C:/Users/facum/OneDrive/Desktop/GWAS/Dataset/BLUPs - Non orthogonalized.xlsx", 1)

#### Bringing the GRM
GRM <- read.delim("C:/Users/facum/OneDrive/Desktop/GWAS/SNP/Kinship.txt", header=FALSE)
lines <- GRM[,1]
GRM <- GRM[-c(1)]
colnames(GRM) <- lines
rownames(GRM) <- lines
GRM <- as.matrix(GRM)
GRM <- make.positive.definite(GRM) # to apply cholesky as suggested by paper implem.

### Defining the trait matrix:
Y <- as.matrix(traits)
n_traits <- ncol(Y)
n_lines <- length(unique(lines))

### Defining CV parameters:
n_partitions <- 20
test_ratio <- 0.2
test_lines <- floor(0.2*n_lines)
train_lines <- n_lines - test_lines

# Initializing each test partition
te_partition_1 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_2 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_3 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_4 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_5 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_6 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_7 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_8 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_9 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_10 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_11 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_12 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_13 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_14 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_15 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_16 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_17 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_18 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_19 <- matrix(NA, nrow=test_lines, n_traits)
te_partition_20 <- matrix(NA, nrow=test_lines, n_traits)

# Initializing each training partition
tr_partition_1 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_2 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_3 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_4 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_5 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_6 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_7 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_8 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_9 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_10 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_11 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_12 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_13 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_14 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_15 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_16 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_17 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_18 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_19 <- matrix(NA, nrow=train_lines, n_traits)
tr_partition_20 <- matrix(NA, nrow=train_lines, n_traits)

Y_test <- list(te_partition_1, te_partition_2, te_partition_3, te_partition_4, te_partition_5,
               te_partition_6, te_partition_7, te_partition_8, te_partition_9,
               te_partition_10, te_partition_11, te_partition_12, te_partition_13,
               te_partition_14, te_partition_15, te_partition_16, te_partition_17,
               te_partition_18, te_partition_19, te_partition_20)

Y_train <- list(tr_partition_1, tr_partition_2, tr_partition_3, tr_partition_4, tr_partition_5,
                tr_partition_6, tr_partition_7, tr_partition_8, tr_partition_9,
                tr_partition_10, tr_partition_11, tr_partition_12, tr_partition_13,
                tr_partition_14, tr_partition_15, tr_partition_16, tr_partition_17,
                tr_partition_18, tr_partition_19, tr_partition_20)

for(j in 1:n_traits){
  
  pheno <- data.frame(GID = lines, Env = '', Response = Y[,j])
  CV_SCHEME <- CV.RandomPart(pheno, NPartitions = n_partitions, PTesting = 0.2, set_seed = 1)

  for (i in 1:n_partitions){
    current_test_partition_indexes <- CV_SCHEME$CrossValidation_list[i]
    current_test_partition_indexes <- as.vector(unlist(current_test_partition_indexes))
    current_test_partition_values <- traits[current_test_partition_indexes, j]
    current_train_partition_values <- traits[-current_test_partition_indexes, j]
    Y_test[[i]][,j] <- current_test_partition_values
    Y_train[[i]][,j] <- current_train_partition_values
  }
}