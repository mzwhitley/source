## reformat GEO dataset, analyzed in Expressionist.
# # 1. primary data is assay id and wide columns of data with sample name as column header
#     melt wide format to long
#     add expt_grp_name from sample key file
#     colnames of final output:
#         expt_grp_name
#         primary_assay_id
#         primary_value
#         replicate_name

# # 2. differential data needs to be log2FC
#     Analyst reports different versions depending on the data started
#       in this case, have  -7 - 7 with full values between -1 and 1 so it looks like a log diff
#     Add dataset_name
#     colnames of final output
#           dataset_name
#           primary_assay_ID
#           reference_id
#           reference_name
#           log2diff
#           pvalue
#           bh



getwd()   
setwd("C:/Users/mwhitley/Box Sync/Functional Genomics/data-AMPK/_20160927_GSE65585/")

in.p <- "PrimaryNormalized_wide.txt"
sample <- "sample.key.txt"
out.p <- "pdata_20160927.txt"

in.d <- "Differential_directedEffect.txt"
char <- c(1:3)
key.d <- c("2016092701")
out.d <- "data_20160927.txt"

library(stringr)
library(reshape2)
write.noRows <- function(input, file1) {write.table(input, file = file1 ,append = FALSE, quote = FALSE, sep = "\t",
                                                    eol = "\n", na = "0", dec = ".", row.names = FALSE,
                                                    col.names = TRUE,
                                                    fileEncoding = "")}
memory.limit(size = 10480)


# read in primary data
#---------------------
data <- read.delim(in.p, sep = "\t")
print(head(data))
# melt
data.m <- melt(data[,c(1:length(colnames(data))-1)])
head(data.m)
# add group_name
sam.key <- read.delim(sample)
data.m <- merge(data.m, sam.key, by.x = "variable", by.y = "Sample", all.x = TRUE)
summary(data.m)
 # sanity check data.m[data.m$Row.Names == "ILMN_2427766",]
 # matched to original GEO (in Expressionist) data
# standard colnames
colnames(data.m) <- c("replicate_name", "primary_assay_id", "primary_value", "expt_grp_name")
# write
write.noRows(data.m, out.p)


# read in differential data
#---------------------

data <- read.delim(in.d, sep = "\t", as.is = char)
print(head(data))
data <- data[,c(1:length(colnames(data))-1)]
# add dataset name
data$dataset_name = key.d[1]
# sanity check: data[data$primary_assay_ID == "ILMN_2783833",]
# print
write.noRows(data, out.d)