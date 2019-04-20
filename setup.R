install.packages("BiocManager", repos = "http://cran.utstat.utoronto.ca/")
install.packages('rjags', repos = "http://cran.utstat.utoronto.ca/") 

library(BiocManager)
install("BSgenome.Hsapiens.UCSC.hg38")
install("MADSEQ") 

library("MADSEQ") 
library("rjags") 
