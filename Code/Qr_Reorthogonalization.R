library(openxlsx)
library(lme4)

rm(list=ls()) 

preprocess_xlsx <- function(filename, sheet){
  file <- read.xlsx(filename, sheet)
  df <- data.frame(file, stringsAsFactors = False)
  headers <- df[,1]
  df <- df[,-1]
  rownames(df) <- headers
  return(df)
}

### Traits
traits <- preprocess_xlsx("/dir/Pheno Traits.xlsx", 1)


for(i in 1:ncol(traits)){
  traits[is.na(traits[,i]), i] <- mean(traits[, i], na.rm = TRUE)
}

### Fixing rank deficiency
traits <- as.matrix(traits)
traits <- traits[,-qr(traits)$pivot[seq_len(qr(traits)$rank)]]
traits <- data.frame(traits)

