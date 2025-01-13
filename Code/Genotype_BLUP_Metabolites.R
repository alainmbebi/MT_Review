library(openxlsx)
library(lme4)

rm(list=ls()) 

preprocess_xlsx <- function(filename, sheet){
  file <- read.xlsx(filename, sheet)
  df <- data.frame(file, stringsAsFactors = False)
  df <- df[,c(-1,-2, -3)]
  return(df)
}

### Metabolic data)
metabolomics <- preprocess_xlsx("/dir/Metabolites.xlsx", 1)

### Replicate data)
replicates <- preprocess_xlsx("/dir/Metabolites.xlsx", 2)

### Nr mets
c_met <-  ncol(replicates)

### Nr genotypes
n_genotypes <- nrow(metabolomics)

### Repeat Index)
I <- c(rep("I", n_genotypes))
II <- c(rep("II", n_genotypes))
Replicate <-  c(I, II)
Replicate <- factor(Replicate, levels = c("I", "II"));

### Genotypes (individuals)
Ind <- c(rep(1, 2*n_genotypes))
for(i in 0:1){
  st <- 1+i*n_genotypes
  fs <- (1+i)*n_genotypes
  Ind[st:fs] <- seq(1, n_genotypes, by=1)
}
Ind <- factor(Ind);


blups = data.frame(matrix(nrow = n_genotypes, ncol = c_met)) 
colnames(blups) <- names(metabolomics)

for(i in  1:c_met){
  
  metabolite <- metabolomics[,i]
  metabolite <- data.frame(metabolite)
  
  replicate <- replicates[,i]
  replicate <- data.frame(replicate)
  
  metabolics <- c(metabolite, replicate)
  metabolics <-  stack(metabolics)
  metabolics <-subset(metabolics, select=values)
  metabolics <- data.frame(metabolics)
  colnames(metabolics) <- "metabolics"
  
  data <- data.frame(cbind(metabolics, Ind, Replicate))
  
  mmodel <- lmer(metabolics ~ (1|Replicate) + (1|Ind)  , data=data, REML=F)  
  
  global_mean <- coef(summary(mmodel))[1]
    
  blup <- ranef(mmodel)
  blup_ind <- data.frame(blup$Ind + global_mean)
  blups[,i] <- blup_ind
}

# Scaling blups min max
blups_scaled <- apply(blups, 2, function(x) {
  (x-min(x))/(max(x)-min(x))
}
)

blups_scaled <- blups_scaled[,colSums(is.na(blups_scaled))<nrow(blups_scaled)]
blups_scaled <- data.frame(blups_scaled)
