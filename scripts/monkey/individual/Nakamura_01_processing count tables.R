library(tidyverse)

### Processing supp tables + processed data from RPM to TPM using reported SC3seq mapped reads ###
ReadsMapped<-as.data.frame(metadata$`SC3seq mapped`)
rownames(ReadsMapped)<-metadata$SampleID
colnames(ReadsMapped)<-'reads'

RPM<-read_tsv('SC3seq_Cy_ProcessedData.txt')
RPM<-as.data.frame(RPM)
gene_ids<-as.character(RPM$macFas5_entrez_id)
gene_names<-RPM$macFas5_gene_symbol
RPM<-RPM[,3:423]
common<-intersect(rownames(ReadsMapped),colnames(RPM))

ReadsMapped<-as.data.frame(ReadsMapped[common,])
rownames(ReadsMapped)<-common
colnames(ReadsMapped)<-'reads'

RPM<-as.matrix(RPM)

Counts<-RPM%*%diag(ReadsMapped$reads)
Counts<-Counts/10e6
colnames(Counts)<-common

### subsetting and re-ordering matrix to match gene lengths file ###

meta<-read.csv('metadata.csv')
meta<-meta[1:390,]
rownames(meta)<-meta$ï..SampleID
emb<-intersect(rownames(meta), colnames(Counts))
Counts<-Counts[,emb]

common_ids<-as.character(intersect(gene_ids,lengths_to_genes$GeneID))
rownames(Counts)<-as.character(gene_ids)
Counts<-Counts[common_ids,]

rownames(lengths_to_genes)<-as.character(lengths_to_genes$GeneID)
lengths_to_genes<-lengths_to_genes[common_ids,]


### TPM normalization ###

tpm <- function(x, gene.length) {
  x <- as.matrix(x)
  len.norm.lib.size <- colSums(x / gene.length)
  return((t(t(x) / len.norm.lib.size) * 1e06) / gene.length)
}

tpm_matrix<-tpm(Counts, lengths_to_genes$length)

rownames(tpm_matrix)<-lengths_to_genes$Gene_name

saveRDS(Counts, 'Counts.rds')
saveRDS(tpm_matrix, 'tpm_matrix.rds')

#### quasi-UMI normalization ####

library(quminorm)


#Now obtaining QUMI counts using custom parameter
#will use shape parameter default of 2.
#Then check if any NA in matrix. If there is, will continue to decrease shape parameter by 0.1 until not.

system.time(qumi_matrix<-quminorm(tpm_matrix, shape=2))
any(is.na(qumi_matrix))

#2 works for no NAs
saveRDS(qumi_matrix, 'qumi_matrix')
