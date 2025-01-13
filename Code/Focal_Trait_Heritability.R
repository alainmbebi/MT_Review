library(heritability)
library(data.table)
library(openxlsx)

setwd("/dir")

# Reading the markers
GRM <- fread("kinship.txt")
GRM <- as.matrix(GRM)
lines <- GRM[, c(1)]
row.names(GRM) <- lines
GRM <- GRM[, -c(1)]
colnames(GRM) <- lines
class(GRM) <- "numeric"

Y <- as.matrix(read.xlsx("Focal.xlsx", sheet=1, rowNames = TRUE))
trait_names <- colnames(Y)

h2_df <- matrix(data=NA, nrow=1, ncol=12)
colnames(h2_df) <- trait_names
for(i in 1:12){
  pheno <- Y[,i]
  pheno <- cbind(lines, pheno)
  colnames(pheno)[1] <- "genotype"
  colnames(pheno)[2] <- trait_names[i]
  pheno <- as.data.frame(pheno)
  class(pheno[,2]) <- "numeric"
  output <- marker_h2(data.vector=pheno[,2], geno.vector=pheno$genotype, K=GRM, h2=TRUE)
  heritability <- output$h2
  h2_df[1, i] <- heritability
}

write.xlsx(as.data.frame(h2_df), "Focal_Traits_h2.xlsx")

