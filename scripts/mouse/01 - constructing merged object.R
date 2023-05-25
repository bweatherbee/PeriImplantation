library(Seurat)
library(tidyverse)
library(biomaRt)
PijuanSala <- readRDS("D:/realigning_for_int/mouse/pijuan_sala/PijuanSala_preorganogenesis.rds")
Mohammed <- readRDS("D:/realigning_for_int/mouse/mohammed/Mohammed_qUMI.RDS")
Deng <- readRDS("D:/realigning_for_int/mouse/Deng/Deng_qUMI.RDS")
Cheng <- readRDS("D:/realigning_for_int/mouse/Cheng/Cheng_qUMI.RDS")

PijuanSala$Paper<-'PijuanSala'
PijuanSala@assays$RNA_qumi<-PijuanSala@assays$UMI_qumi
PijuanSala@assays$UMI_qumi<-NULL
PijuanSala@active.assay<-'RNA_qumi'
PijuanSala$Stage<-PijuanSala$stage
PijuanSala$stage<-NULL
merged<-merge(x=Deng, y=c(Mohammed, PijuanSala, Cheng))

Idents(merged)<-merged$Age
merged<-subset(merged, idents=c('D03.5', 'D04', 'D04.5', 'D05.25', 'D05.5', 'D06.25', 'D06.5', 'D06.75', 'D07.0', 'D07.25', 'D07.5'))


merged_list<-SplitObject(merged, split.by="Paper")
merged_list<-lapply(X=merged_list, FUN=SCTransform, assay='RNA_qumi')

features<-SelectIntegrationFeatures(object.list=merged_list, nfeatures=5000)
merged_list<-PrepSCTIntegration(object.list=merged_list, anchor.features = features)
anchors<-FindIntegrationAnchors(object.list=merged_list,normalization.method = 'SCT', anchor.features=features)
merged<-IntegrateData(anchorset=anchors, normalization.method = 'SCT')


### importing methods from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/ ###


# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  
  return(humanx)
}

m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
merged<-CellCycleScoring(merged, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = F)

DefaultAssay(merged)<-"RNA_qumi"
merged<-ScaleData(merged, vars.to.regress=c("S.Score", "G2M.Score", "percent.mt"),features=rownames(merged))

DefaultAssay(merged)<-"integrated"
merged<-ScaleData(merged, vars.to.regress=c("S.Score", "G2M.Score", "percent.mt"),features=rownames(merged))

merged<-RunPCA(merged, features=rownames(merged))
ElbowPlot(merged, ndims = 50)
merged<-RunUMAP(merged, dims=1:10)
merged<-FindNeighbors(merged, reduction="pca", dims=1:10)
merged<-FindClusters(merged, resolution = 0.02)
DimPlot(merged)
DimPlot(merged, group.by='Paper')
DimPlot(merged, group.by='Age')
DimPlot(merged, group.by='Stage')

DefaultAssay(merged)<-'RNA_qumi'
FeaturePlot(merged, features=c('Pou5f1', 'Nanog', 'Sox2',
                               'Gata6', 'Sox17', 'T',
                               'Gata3', 'Cdx2', 'Eomes'), order=T, min.cutoff = 'q20', max.cutoff = 'q98')

new.cluster.ids<-c('EPI', 'EXE', 'VE')
names(new.cluster.ids)<-levels(merged)
merged<-RenameIdents(merged, new.cluster.ids)
DimPlot(merged)


merged$new_cell_type<-Idents(merged)

metadata<-merged@meta.data


df<-metadata$Age
df<-as.data.frame(df)
df$new_cell_type<-metadata$new_cell_type
df$Stage<-metadata$Stage
rownames(df)<-rownames(metadata)
colnames(df)<-c("Age", "new_cell_type", "Stage")
df$new_cell_type<-as.character(df$new_cell_type)

df <- within(df, {
  f <- Age == 'D06.5'
  Age[f] <- 'D06.5'
  Stage[f] <- 'late_post'
}) 

df <- within(df, {
  f <- Age == 'D03.5'
  Age[f] <- 'D03.5'
  Stage[f] <- 'blastocyst'
}) 

df <- within(df, {
  f <- Age == 'D04'
  Age[f] <- 'D014'
  Stage[f] <- 'blastocyst'
}) 

df <- within(df, {
  f <- Age == 'D04.5'
  Age[f] <- 'D04.5'
  Stage[f] <- 'blastocyst'
})

merged$new_cell_type<-df$new_cell_type
merged$Stage<-df$Stage

DimPlot(merged, group.by = 'new_cell_type')
DimPlot(merged, group.by='Stage')
DimPlot(merged,group.by='Age')
Idents(merged)<-merged$new_cell_type
DefaultAssay(merged)<-'RNA_qumi'
merged<-SCTransform(merged, assay='RNA_qumi')
saveRDS(merged, 'merged_mouse.rds')

df<-metadata$age
df<-as.data.frame(df)
df$Stage<-metadata$Stage
rownames(df)<-rownames(metadata)
colnames(df)<-c("Age","Stage")

df <- within(df, {
  f <- Stage == 'early_blast'
  Stage[f] <- 'blastocyst'
}) 

df <- within(df, {
  f <- Stage == 'late_blast'
  Stage[f] <- 'blastocyst'
}) 


my_levels<-c("blastocyst", "early_post", "late_post", "gast")
merged$Stage<- factor(x = merged$Stage, levels = my_levels)
f
my_levels<-c("EPI", "VE", "EXE")
merged$new_cell_type<-factor(x = merged$new_cell_type, levels = my_levels)

saveRDS(merged, 'merged_mouse.rds')


#### converting to human gene names ####

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  
  return(genesV2)
}


DefaultAssay(merged)<-'RNA_qumi'
RNA_qumi_genes<-rownames(merged)
DefaultAssay(merged)<-'SCT'
SCT_genes<-rownames(merged)
DefaultAssay(merged)<-'integrated'
integrated_genes<-rownames(merged)


key_RNA_qumi<-convertMouseGeneList(RNA_qumi_genes)
key_SCT<-convertMouseGeneList(SCT_genes)
key_integrated<-convertMouseGeneList(integrated_genes)

DefaultAssay(merged)<-'RNA_qumi'
merged_hgnc<-merged[key_RNA_qumi$MGI.symbol,]
rownames(key_RNA_qumi)<-make.unique(key_RNA_qumi$MGI.symbol)
key_RNA_qumi<-key_RNA_qumi[rownames(merged_hgnc),]
merged_hgnc@assays[["RNA_qumi"]]@counts@Dimnames[[1]]<-key_RNA_qumi$HGNC.symbol
merged_hgnc@assays[["RNA_qumi"]]@data@Dimnames[[1]]<-key_RNA_qumi$HGNC.symbol

DefaultAssay(merged_hgnc)<-'SCT'
rownames(key_SCT)<-make.unique(key_SCT$MGI.symbol)
key_SCT<-key_SCT[rownames(merged_hgnc),]
merged_hgnc@assays[["SCT"]]@counts@Dimnames[[1]]<-key_SCT$HGNC.symbol
merged_hgnc@assays[["SCT"]]@data@Dimnames[[1]]<-key_SCT$HGNC.symbol

DefaultAssay(merged_hgnc)<-'integrated'
rownames(key_integrated)<-make.unique(key_integrated$MGI.symbol)
key_integrated<-key_integrated[rownames(merged_hgnc),]
merged_hgnc@assays[["integrated"]]@data@Dimnames[[1]]<-key_integrated$HGNC.symbol

saveRDS(key_RNA_qumi, 'key_mouse_to_human.rds')
saveRDS(merged_hgnc, 'merged_mouse_hgncsymbols.rds')

### DEG analysis ###

DefaultAssay(merged)<-'SCT'
Idents(merged)<-merged$new_cell_type
merged<-PrepSCTFindMarkers(merged)

lineage_DEGs<-FindAllMarkers(merged, assay='SCT', logfc.threshold=0, test.use='roc', only.pos=T)

epi<-subset(merged, idents='EPI')
ve<-subset(merged, idents='VE')
exe<-subset(merged, idents='EXE')

Idents(epi)<-epi$Stage
Idents(ve)<-ve$Stage
Idents(exe)<-exe$Stage

epi<-PrepSCTFindMarkers(epi)
ve<-PrepSCTFindMarkers(ve)
exe<-PrepSCTFindMarkers(exe)

epi_stage<-FindAllMarkers(epi, assay='SCT', logfc.threshold=0, test.use='roc', only.pos=T)
ve_stage<-FindAllMarkers(ve, assay='SCT', logfc.threshold=0, test.use='roc', only.pos=T)
exe_stage<-FindAllMarkers(exe, assay='SCT', logfc.threshold=0, test.use='roc', only.pos=T)

write.csv(lineage_DEGs, 'tables/lineage_DEGs.csv')
write.csv(epi_stage, 'tables/epi_stage_DEGs.csv')
write.csv(ve_stage, 'tables/ve_stage_DEGs.csv')
write.csv(exe_stage, 'tables/exe_stage_DEGs.csv')

average_expression<-as.data.frame(AverageExpression(merged, assays='SCT', add.ident='Stage'))
write.csv(average_expression, 'tables/mouse_avgSCTExpression_stage_celltype.csv')
