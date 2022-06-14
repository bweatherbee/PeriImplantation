library(tidyverse)
library(Seurat)

Xiang_metadata <- read_csv("C:/Users/baile/Desktop/realigning_for_int/Xiang/Xiang_metadata.csv")
Xiang_metadata<-as.data.frame(Xiang_metadata)
rownames(Xiang_metadata)<-Xiang_metadata$Cell_id
Xiang_metadata<-Xiang_metadata[,2:5]
matrix<-readRDS('matrix.RDS')
tpm_matrix<-readRDS('tpm_matrix.RDS')
qumi_matrix<-readRDS('qumi_matrix.RDS')

Xiang_metadata<-Xiang_metadata[colnames(qumi_matrix),]



Xiang<-CreateSeuratObject(counts=qumi_matrix, project="Xiang",
                          assay="RNA_qumi", meta.data = Xiang_metadata)
Xiang@assays[["RNA_TPM"]]<-CreateAssayObject(counts = tpm_matrix)
Xiang@assays[["RNA_estCounts"]]<-CreateAssayObject(counts=matrix)

Xiang_TPM<-CreateSeuratObject(counts=tpm_matrix, project="Xiang",
                          assay="RNA_tpm", meta.data = Xiang_metadata)
Xiang_TPM[['percent.mt']]<-PercentageFeatureSet(Xiang_TPM, pattern="^MT")

VlnPlot(Xiang_TPM, features=c('nFeature_RNA_tpm', 'nCount_RNA_tpm', 'percent.mt'))


DefaultAssay(Xiang)<-'RNA_qumi'

Xiang[['percent.mt']]<-PercentageFeatureSet(Xiang, pattern="^MT")

VlnPlot(Xiang, features=c('nFeature_RNA_qumi', 'nCount_RNA_qumi', 'percent.mt'))

Xiang<-NormalizeData(Xiang, normalization.method = "RC", scale.factor = 10e6)
Xiang<-FindVariableFeatures(Xiang)
all.genes<-rownames(Xiang)
Xiang<-ScaleData(Xiang, features=all.genes)

Xiang<-RunPCA(Xiang, features=VariableFeatures(object=Xiang))
DimHeatmap(Xiang, dims=1:15, cells=500, balanced=TRUE)

ElbowPlot(Xiang, ndims=50)

Xiang<-FindNeighbors(Xiang, dims=1:5)
Xiang<-FindClusters(Xiang, resolution=0.2)

Xiang<-RunUMAP(Xiang, dims=1:20)
DimPlot(Xiang, reduction="umap")
DimPlot(Xiang, group.by = "Age")

FeaturePlot(Xiang, features=c("NANOG", "SOX2", "POU5F1","TBXT",
                              "GATA4", "SOX17", "GATA6", "PDGFRA",
                              "GATA3", "GATA2", "SDC1", "HLA-G"),
            order=TRUE, slot='counts')


new.cluster.ids<-c("TrB", "EPI", "STB", "PSA-EPI", "HYPO", "EVT")
names(new.cluster.ids)<-levels(Xiang)
Xiang<-RenameIdents(Xiang, new.cluster.ids)
DimPlot(Xiang)
Xiang$cell_type<-Idents(Xiang)
saveRDS(Xiang, "Xiang_qUMI.RDS")
