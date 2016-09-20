setwd("C:/Users/mwhitley/Box Sync")

# source("https://bioconductor.org/biocLite.R")
# biocLite("ArrayExpress")
# biocLite("affy")


library(ArrayExpress)
library(affy)


load("C:/Users/mwhitley/Downloads/E-MEXP-3291.eSet.R")
norm <- rma(study)
matrix <- exprs(norm)
write.table(pData(norm), "E-MEXP-3291.pData.txt", sep = "\t", row.names = TRUE, col.names = NA) 
write.table(matrix,"E-MEXP-3291.exprs.txt", sep = "\t", row.names = TRUE, col.names = NA)
