###Question 1
#Inside of the supercomputer
#ls

###Question 2
#inside of R
#Install ShortRead package to load in fastq files
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead")
#Load package
library(ShortRead)
sequences <- readFastq(dirPath = "fastq/SRR5324768_pass_1.fastq")
sequences
#length: 250803 reads; width: 101

