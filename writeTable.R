
write <- function(input, file1) { 
  
  write.table(input, file = file1 ,append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "0", dec = ".", row.names = TRUE,
            col.names = NA,
            fileEncoding = "")
}
