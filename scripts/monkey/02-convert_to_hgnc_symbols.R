library(tidyverse)
t2g<-read_tsv('Cyn5.0_t2g.txt', col_names = F)
t2g<-t2g[,2:3]
colnames(t2g)<-c('geneID', 'gene_name')
t2g<-unique(t2g)

DefaultAssay(merged)<-'SCT'
rows_SCT<-rownames(merged)
DefaultAssay(merged)<-'RNA'
rows_RNA<-rownames(merged)

t2g<-as.data.frame(t2g)
rownames(t2g)<-make.unique(t2g$gene_name)
map_RNA<-t2g[rows_RNA,]
map_SCT<-t2g[rows_SCT,]

human_to_macFas <- read_excel("human_cynmonkey_genemaps.xlsx")

human_to_macFas<-as.data.frame(human_to_macFas)
rownames(human_to_macFas)<-as.character(human_to_macFas$macFas5_entrez_id)

rownames(map_RNA)<-map_RNA$geneID
rownames(map_SCT)<-map_SCT$geneID

common_RNA<-intersect(rownames(human_to_macFas), rownames(map_RNA))
common_SCT<-intersect(rownames(human_to_macFas), rownames(map_SCT))

DefaultAssay(merged)<-'SCT'
merged@assays[["SCT"]]@counts@Dimnames[[1]]<-rownames(map_SCT)
merged@assays[["SCT"]]@data@Dimnames[[1]]<-rownames(map_SCT)

DefaultAssay(merged)<-'RNA'
merged@assays[["RNA"]]@counts@Dimnames[[1]]<-rownames(map_RNA)
merged@assays[["RNA"]]@data@Dimnames[[1]]<-rownames(map_RNA)

key_RNA<-human_to_macFas[common_RNA,]
key_SCT<-human_to_macFas[common_SCT,]

DefaultAssay(merged)<-'RNA'

subset_1<-subset(merged, features = common_RNA)
subset_1@assays[["RNA"]]@counts@Dimnames[[1]]<-key_RNA$hg19_gene_symbol
subset_1@assays[["RNA"]]@data@Dimnames[[1]]<-key_RNA$hg19_gene_symbol
subset_1@assays[["SCT"]]@counts@Dimnames[[1]]<-key_SCT$hg19_gene_symbol
subset_1@assays[["SCT"]]@data@Dimnames[[1]]<-key_SCT$hg19_gene_symbol

merged_hgnc<-subset_1
saveRDS(merged_hgnc, 'merged_hgnc19_symbols.RDS')
