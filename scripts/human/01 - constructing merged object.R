library(Seurat)
library(data.table)
library(tidyverse)

Petropoulos <- readRDS("D://realigning_for_int/Petropoulos/Petropoulos_qUMI.RDS")
MZG <- readRDS("D://realigning_for_int/MZG/MZG.RDS")
MZG$percent.mt<-MZG$percent.mito
MZG$species<-"human"
MZG$paper<-"Mole"
DefaultAssay(MZG)<-"RNA"
MZG$percent.mt<-MZG$percent.mito
MZG$percent.mito<-NULL
MZG$nCount_SCENIC<-NULL
MZG$nFeature_SCENIC<-NULL
MZG@assays$SCENIC<-NULL
MZG@assays$dorothea<-NULL
MZG<-NormalizeData(MZG, normalization.method="RC", scale.factor = 10e6)

#Yan <- readRDS("D://realigning_for_int/Yan/Yan_qUMI.RDS") # excluding Yan here because blastocyst cells are <50 total - makes integration difficult
Blakely <- readRDS("D://realigning_for_int/Blakely/Blakely_qUMI.RDS")
Zhou <- readRDS("D://realigning_for_int/MZG_Xiang_Zhou_VERSION2/datasets/Zhou_UMIcounts.rds")
Zhou$Age<-Zhou$timepoint
Zhou$timepoint<-NULL
Zhou<-NormalizeData(Zhou, normalization.method = "RC", scale.factor = 10e6)

Xiang <- readRDS("D://realigning_for_int/Xiang/gene_matrix_estimated_counts_output/August2021/Xiang_qUMI.RDS")


Blakely_matrix<-Blakely@assays[["RNA_qumi"]]@counts
Blakely_meta<-Blakely@meta.data
MZG_matrix<-MZG@assays[["RNA"]]@counts
MZG_meta<-MZG@meta.data
Petropoulos_matrix<-Petropoulos@assays[["RNA_qumi"]]@counts
Petropoulos_meta<-Petropoulos@meta.data
Xiang_matrix<-Xiang@assays[["RNA_qumi"]]@counts
Xiang_meta<-Xiang@meta.data
Zhou_matrix<-Zhou@assays[["RNA"]]@counts
Zhou_meta<-Zhou@meta.data

common.features<-intersect(rownames(Blakely_matrix), rownames(MZG_matrix))
common.features<-intersect(rownames(Zhou_matrix), common.features)
Blakely_matrix<-Blakely_matrix[common.features,]
MZG_matrix<-MZG_matrix[common.features,]
Petropoulos_matrix<-Petropoulos_matrix[common.features,]
Xiang_matrix<-Xiang_matrix[common.features,]
Zhou_matrix<-Zhou_matrix[common.features,]

merged_matrix<-merge(x=Blakely_matrix, y=c(MZG_matrix, Petropoulos_matrix, Xiang_matrix, Zhou_matrix),by=0)
rownames(merged_matrix)<-merged_matrix$Row.names
merged_matrix<-merged_matrix[,2:10135]

common.meta<-intersect(colnames(Blakely_meta), colnames(MZG_meta))
Blakely_meta<-Blakely_meta[,common.meta]
MZG_meta<-MZG_meta[,common.meta]
Petropoulos_meta$cell_type<-"TBC"
Petropoulos_meta<-Petropoulos_meta[,common.meta]
Xiang_meta<-Xiang_meta[,common.meta]
Zhou_meta$percent.mt<-"TBC"
Zhou_meta$seurat_clusters<-"TBC"
Zhou_meta<-Zhou_meta[,common.meta]

l=list(Blakely_meta, MZG_meta, Petropoulos_meta, Xiang_meta, Zhou_meta)
merged_meta<-rbindlist(l, use.names=TRUE)
rownames(merged_meta)<-colnames(merged_matrix)

merged<-CreateSeuratObject(counts=merged_matrix, project="Human_merged", assay="UMI_qumi",
                           meta.data=merged_meta)

merged[['percent.mt']]<-PercentageFeatureSet(merged, pattern="^MT")
VlnPlot(merged, features=c('nFeature_UMI_qumi', 'nCount_UMI_qumi', 'percent.mt'))

merged<-NormalizeData(merged, normalization.method = "RC", scale.factor = 10e6)

merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D1", "D01"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D2", "D02"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D3", "D03"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D4", "D04"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D5", "D05"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D6", "D06"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D7", "D07"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D010", "D10"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D011", "D11"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D012", "D12"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D014", "D14"))

#Realized there is a repeat from the end of Petropoulos, removing
cells<-colnames(merged)
removeWords <- function(str, stopwords) {
  x <- unlist(strsplit(str, " "))
  paste(x[!x %in% stopwords])
}
cells<-removeWords(cells, 'E3.53.3438.1')

merged<-subset(merged, cells=cells)

### For signaling paper, only interested in blastocyst stage on so subsetting for those ages ###
Idents(merged)<-merged$Age

merged<-subset(merged, idents=c('D05' ,'D06', 'D07', 'D08', 'D09', 'D10', 'D11', 'D12', 'D14'))



#### Trying cell cycle and mitochondrial content regressions + Merging Using SCTransform ###

DefaultAssay(merged)<-"UMI_qumi"
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
merged<-CellCycleScoring(merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(merged, features=c("PCNA", "TOP2A", "MCM6", "MKI67", ncol=2))
VlnPlot(merged, features="percent.mt")
merged<-ScaleData(merged, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
saveRDS(merged, 'combined_scaled.rds')

#merged<-RenameAssays(merged, 'UMI_qumi'='RNA')
merged_list<-SplitObject(merged, split.by="paper")
merged_list<-lapply(X=merged_list, FUN=SCTransform, assay='UMI_qumi')

features<-SelectIntegrationFeatures(object.list=merged_list, nfeatures=3000)
merged_list<-PrepSCTIntegration(object.list=merged_list, anchor.features=features)
anchors<-FindIntegrationAnchors(object.list=merged_list, normalization.method='SCT', anchor.features=features)
merged<-IntegrateData(anchorset=anchors,k.weight=80)
saveRDS(merged, 'merged_sctint.RDS')

DefaultAssay(merged)<-"integrated"
merged<-ScaleData(merged, vars.to.regress=c("S.Score", "G2M.Score", "percent.mt"),features=rownames(merged))
merged<-RunPCA(merged, features=rownames(merged))
ElbowPlot(merged, ndims = 50)
merged<-RunUMAP(merged, dims=1:30)
merged<-FindNeighbors(merged, reduction="pca", dims=1:10)
merged<-FindClusters(merged, resolution = 0.2)
DimPlot(merged)

DimPlot(merged, group.by = 'paper')
DimPlot(merged, group.by='Age')
DimPlot(merged)

DefaultAssay(merged)<-"UMI_qumi"
FeaturePlot(merged, features=c("SOX2", "NANOG", "POU5F1", "KLF17"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("SOX17", "GATA4", "GATA6", "PDGFRA"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("GATA3", "GATA2","TFAP2C", "CDX2"),order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("SDC1", "CGB1", "CGA", "ITGB3"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')


new.cluster.ids<-c("TE", "STB", "EVT", "TE", "EPI", "TE", "HYPO", "TE", "TE")
names(new.cluster.ids)<-levels(merged)
merged<-RenameIdents(merged, new.cluster.ids)
DimPlot(merged)

merged$new_cell_type<-Idents(merged)

metadata<-merged@meta.data
df<-metadata$Age
df<-as.data.frame(df)
df$new_cell_type<-metadata$new_cell_type
df$stage<-metadata$stage
rownames(df)<-rownames(metadata)
colnames(df)<-c("Age", "new_cell_type", "stage")
df$new_cell_type<-as.character(df$new_cell_type)

df <- within(df, {
  f <- stage == 'blastocyst' & new_cell_type=='STB'
  stage[f] <- 'blastocyst'
  new_cell_type[f] <- 'TE'
}) 

df <- within(df, {
  f <- stage == 'blastocyst' & new_cell_type=='EVT'
  stage[f] <- 'blastocyst'
  new_cell_type[f] <- 'TE'
}) 

merged$new_cell_type<-df$new_cell_type
Idents(merged)<-merged$new_cell_type

levels(merged$new_cell_type)
my_levels<-c("EPI", "HYPO", "TE", "STB", "EVT")
merged$new_cell_type <- factor(x = merged$new_cell_type, levels = my_levels)

my_levels<-c("blastocyst", "peri", "early post", "late post")
merged$stage<- factor(x = merged$stage, levels = my_levels)

merged<-saveRDS(merged, 'merged_integrated.RDS')

### Differential expression analyses ###

DefaultAssay(merged)<-'SCT'
Idents(merged)
merged_3lin<-subset(merged, idents=c('EPI', 'HYPO', 'TE'))
merged_3lin<-PrepSCTFindMarkers(merged_3lin)
lineage_DEGs<-FindAllMarkers(merged_3lin, assay='SCT', logfc.threshold = 0, test.use='roc', min.pct=0.05, only.pos = T)

merged_TrBsubset<-subset(merged, idents=c('TE', 'STB', 'EVT'))
merged_TrBsubset<-PrepSCTFindMarkers(merged_TrBsubset)
TrB_DEGs<-FindAllMarkers(merged_TrBsubset, assay='SCT', logfc.threshold = 0, test.use='roc', min.pct=0.05, only.pos = T)

epi<-subset(merged, idents='EPI')
hypo<-subset(merged, idents='HYPO')
te<-subset(merged, idents='TE')

Idents(epi)<-epi$stage
Idents(hypo)<-hypo$stage
Idents(te)<-te$stage

epi<-PrepSCTFindMarkers(epi)
hypo<-PrepSCTFindMarkers(hypo)
te<-PrepSCTFindMarkers(te)

epi_stages<-FindAllMarkers(epi, assay='SCT', logfc.threshold = 0, test.use='roc', min.pct=0.05, only.pos = T)
hypo_stages<-FindAllMarkers(hypo, assay='SCT', logfc.threshold = 0, test.use='roc', min.pct=0.05, only.pos = T)
te_stages<-FindAllMarkers(te, assay='SCT', logfc.threshold = 0, test.use='roc', min.pct=0.05, only.pos = T)

write.csv(lineage_DEGs, 'tables/epi_hypo_te_DEGs.csv')
write.csv(TrB_DEGs, 'tables/TrB_terminaldiffs_DEGs.csv')
write.csv(epi_stages, 'tables/Epi_stages_DEGs.csv')
write.csv(hypo_stages, 'tables/Hypo_stages_DEGs.csv')
write.csv(te_stages, 'tables/TE_stages_DEGs.csv')

average_expression<-as.data.frame(AverageExpression(merged, assays='SCT', add.ident='stage'))
write.csv(average_expression, 'tables/human_avgSCTExpression_stage_celltype.csv')
