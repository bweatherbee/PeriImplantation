library(tidyverse)


#### Generating matrix from kallisto/kb-python csvs + tpm conversion ####

X<-read.csv("X.csv", header=FALSE)
var<-read.csv("var.csv")
obs<-read.csv("obs.csv")

X<-as.data.frame(X)

rownames(X)<-make.unique(obs$cell_id)
colnames(X)<-var$gene_name

matrix<-t(X)

#tpm matrix function taken from ribiosNGS package
tpm <- function(x, gene.length) {
  x <- as.matrix(x)
  len.norm.lib.size <- colSums(x / gene.length)
  return((t(t(x) / len.norm.lib.size) * 1e06) / gene.length)
}

tpm_matrix<-tpm(matrix, var$length)

saveRDS(matrix, "matrix.rds")
saveRDS(tpm_matrix, "tpm_matrix.rds")


#### quasi-UMI normalization ####

library(quminorm)


#Now obtaining QUMI counts using custom parameter
#Using 2.945 results in mostly NA values. So instead, will use shape parameter default of 2.
#Then check if any NA in matrix. If there is, will continue to decrease shape parameter by 0.1 until not.

system.time(qumi_matrix<-quminorm(tpm_matrix, shape=1.7))
any(is.na(qumi_matrix))

#1.7 accomplished no NA for Blakely Matrix

hist(log1p(rowSums(qumi_matrix)), main="Blakely log(1+QUMI counts)")
hist(log1p(rowSums(matrix)), main="Blakely log(1+est counts)")

saveRDS(qumi_matrix, "qumi_matrix.rds")
