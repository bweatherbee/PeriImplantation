library(Seurat)

merged<-merge(Nakamura, c(Ma, Yang))
merged.list<-SplitObject(merged, split.by = "paper")
merged.list<-lapply(merged.list, FUN=function(x) {
  x<-NormalizeData(x, normalization.method = 'RC', scale.factor = 10e6)
  x<-FindVariableFeatures(x, selection.method = 'vst', nfeatures=2000)
})
features<-SelectIntegrationFeatures(merged.list)
anchors<-FindIntegrationAnchors(merged.list, anchor.features = features)

merged<-IntegrateData(anchors)

DefaultAssay(merged)<-'integrated'
merged<-ScaleData(merged)
merged<-RunPCA(merged, npcs=30)
ElbowPlot(merged)
merged<-RunUMAP(merged, reduction='pca', dims=1:10)
merged<-FindNeighbors(merged, reduction='pca', dims=1:10)
merged<-FindClusters(merged, resolution=0.02)


DimPlot(merged)
DimPlot(merged, group.by = 'Age')
N<-subset(merged, (paper=='Nakamura'))
DimPlot(merged, group.by='paper')
DimPlot(N, group.by='cell_type')

merged<-RenameIdents(merged, `0`="TE", `1`="EPI", `2`="EXMC", `3`="HYPO")
DimPlot(merged)
merged$new_cell_type<-Idents(merged)

DefaultAssay(merged)<-"RNA"
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
merged<-CellCycleScoring(merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(merged, features=c("PCNA", "TOP2A", "MCM6", "MKI67", ncol=2))
VlnPlot(merged, features="percent.mt")
merged<-ScaleData(merged, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))

merged<-SCTransform(merged, method="glmGamPoi", vars.to.regress = c('percent.mt', 'S.Score', 'G2M.Score'))

saveRDS(merged, 'merged_cynmonkey.rds')

Idents(merged)<-merged$new_cell_type
DefaultAssay(merged)<-'SCT'
FeaturePlot(merged, features=c("ISL1", "GABRP", "VTCN1", "DLK1"), min.cutoff = 'q1', max.cutoff = 'q99', order=T)
FeaturePlot(merged, features=c("EOMES","T", "MIXL1", "FOXH1"), min.cutoff = 'q1', max.cutoff = 'q99', order=T)
FeaturePlot(merged, features=c("GATA6", "GATA4", "SOX17", "COL6A1"),min.cutoff = 'q1', max.cutoff = 'q99', order=T)
FeaturePlot(merged, features=c("KLF17", "POU5F1", "SOX2", "NANOG"),min.cutoff = 'q1', max.cutoff = 'q99', order=T)
FeaturePlot(merged, features=c("KRT7", "KRT8", "KRT18", "KRT19"),min.cutoff = 'q20', order=T)
FeaturePlot(merged, features=c("GATA3", "TFAP2C", "TFAP2A", "GATA2"),min.cutoff = 'q1', max.cutoff = 'q99', order=T)
FeaturePlot(merged, features=c("CER1", "LEFTY1", "NOG", "DKK1"),min.cutoff = 'q10', order=T)

DimPlot(merged, group.by='stage')

merged<-PrepSCTFindMarkers(merged)
lineagedegs<-FindAllMarkers(merged, assay='SCT', logfc.threshold = 0, test.use='roc', min.pct=.05, only.pos=T)

epi<-subset(merged, idents='EPI')
hypo<-subset(merged, idents='HYPO')
te<-subset(merged, idents='TE')

Idents(epi)<-epi$stage
Idents(hypo)<-hypo$stage
Idents(te)<-te$stage

epi<-PrepSCTFindMarkers(epi)
hypo<-PrepSCTFindMarkers(hypo)
te<-PrepSCTFindMarkers(te)

epi_stage<-FindAllMarkers(epi, assay='SCT', logfc.threshold = 0, test.use='roc', min.pct=.05, only.pos=T)
hypo_stage<-FindAllMarkers(hypo, assay='SCT', logfc.threshold = 0, test.use='roc', min.pct=.05, only.pos=T)
te_stage<-FindAllMarkers(te, assay='SCT', logfc.threshold = 0, test.use='roc', min.pct=.05, only.pos=T)

write.csv(lineagedegs, 'tables/all_lineages_DEGs.csv')
write.csv(epi_stage, 'tables/epi_stage_DEGs.csv')
write.csv(hypo_stage, 'tables/hypo_stage_DEGs.csv')
write.csv(te_stage, 'tables/te_stage_DEGs.csv')

average_expression<-as.data.frame(AverageExpression(merged, assays='SCT', add.ident='stage'))
write.csv(average_expression, 'tables/monkey_avgSCTExpression_stage_celltype.csv')
