library(tidyverse)
library(Seurat)

Blakely_metadata <- read_excel("C:/Users/baile/Desktop/realigning_for_int/Blakely/Blakely_metadata.xlsx")
Blakely_metadata<-as.data.frame(Blakely_metadata)
rownames(Blakely_metadata)<-make.unique(Blakely_metadata$cell_name)
Blakely_metadata<-Blakely_metadata[,2:6]
matrix<-readRDS('matrix.RDS')
tpm_matrix<-readRDS('tpm_matrix.RDS')
qumi_matrix<-readRDS('qumi_matrix.RDS')


Blakely<-CreateSeuratObject(counts=qumi_matrix, project="Blakely",
                          assay="RNA_qumi", meta.data = Blakely_metadata)
Blakely@assays[["RNA_TPM"]]<-CreateAssayObject(counts = tpm_matrix)
Blakely@assays[["RNA_estCounts"]]<-CreateAssayObject(counts=matrix)

Blakely_TPM<-CreateSeuratObject(counts=tpm_matrix, project="Blakely",
                              assay="RNA_tpm", meta.data = Blakely_metadata)
Blakely_TPM[['percent.mt']]<-PercentageFeatureSet(Blakely_TPM, pattern="^MT")

VlnPlot(Blakely_TPM, features=c('nFeature_RNA_tpm', 'nCount_RNA_tpm', 'percent.mt'))

DefaultAssay(Blakely)<-'RNA_qumi'

Blakely[['percent.mt']]<-PercentageFeatureSet(Blakely, pattern="^MT")

VlnPlot(Blakely, features=c('nFeature_RNA_qumi', 'nCount_RNA_qumi', 'percent.mt'))

Blakely<-NormalizeData(Blakely, normalization.method = "RC", scale.factor = 10e6)
Blakely<-FindVariableFeatures(Blakely)
all.genes<-rownames(Blakely)
Blakely<-ScaleData(Blakely, features=all.genes)

Blakely<-RunPCA(Blakely, features=VariableFeatures(object=Blakely))
DimHeatmap(Blakely, dims=1:15, cells=500, balanced=TRUE)

ElbowPlot(Blakely, ndims=50)

Blakely<-FindNeighbors(Blakely, dims=1:50)
Blakely<-FindClusters(Blakely, resolution=1)

Blakely<-RunUMAP(Blakely, dims=1:50)
DimPlot(Blakely, reduction="umap")
DimPlot(Blakely, group.by = "Age")
DimPlot(Blakely, group.by="orig_cell_type")

FeaturePlot(Blakely, features=c("NANOG", "SOX2", "POU5F1","TBXT",
                              "GATA4", "SOX17", "GATA6", "PDGFRA",
                              "GATA3", "GATA2", "ZSCAN10", ""),
            order=TRUE, slot='counts')

Idents(Blakely)<-Blakely$orig_cell_type
Blakely$cell_type<-Idents(Blakely)
saveRDS(Blakely, "Blakely_qUMI.RDS")
