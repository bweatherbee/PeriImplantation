library(tidyverse)


#### Generating matrix from kallisto/kb-python csvs + tpm conversion ####

X<-read.csv("X.csv", header=FALSE)
var<-read.csv("var.csv")
obs<-read.csv("obs.csv")

X<-as.data.frame(X)

rownames(X)<-obs$cell_id
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

remotes::install_github("willtownes/quminorm")
library(quminorm)

#testing poisson fit on MZG UMI count data to find shape parameter
MZG <- readRDS("C:/Users/baile/Desktop/realigning_for_int/MZG/MZG.RDS")
MZG<-subset(MZG, idents=c("Epiblast", "Hypoblast", "Syncytiotrophoblasts"))
MZG_matrix<-MZG@assays[["RNA"]]@counts
MZG_matrix<-as.matrix(MZG_matrix)

m<-rowSums(MZG_matrix)
keep<-which(colSums(MZG_matrix)>5000)

fit<-poilog_mle(MZG_matrix[,sample(keep,)])
summary(fit$sig)

#Shape parameter median is 2.945 and mean is 3.037


hist(log1p(m),main="MZG log(1+UMI counts)")


#Now obtaining QUMI counts using custom parameter
#Using 2.945 results in mostly NA values. So instead, will use shape parameter default of 2.
#Will then check if any NA in matrix. then will decrease by 0.1 until none.

system.time(qumi_matrix<-quminorm(tpm_matrix, shape=1.6))
any(is.na(qumi_matrix))
#no NAs in Petropolous dataset with shape parameter of 1.6.

hist(log1p(rowSums(qumi_matrix)), main="Petropoulos log(1+QUMI counts)")
hist(log1p(rowSums(matrix)), main="Yan log(1+est counts)")

saveRDS(qumi_matrix, "qumi_matrix.rds")
