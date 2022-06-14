library(tidyverse)
library(Seurat)

Petropoulos_metadata <- read_excel("Petropoulos_metadata.xlsx")
Petropoulos_metadata<-as.data.frame(Petropoulos_metadata)
rownames(Petropoulos_metadata)<-make.unique(Petropoulos_metadata$cell_name)
Petropoulos_metadata<-Petropoulos_metadata[,2:5]
matrix<-readRDS('matrix.RDS')
tpm_matrix<-readRDS('tpm_matrix.RDS')
qumi_matrix<-readRDS('qumi_matrix.RDS')


Petropoulos<-CreateSeuratObject(counts=qumi_matrix, project="Petropoulos",
                          assay="RNA_qumi", meta.data = Petropoulos_metadata)
Petropoulos@assays[["RNA_TPM"]]<-CreateAssayObject(counts = tpm_matrix)
Petropoulos@assays[["RNA_estCounts"]]<-CreateAssayObject(counts=matrix)

Petropoulos_TPM<-CreateSeuratObject(counts=tpm_matrix, project="Petropoulos",
                              assay="RNA_tpm", meta.data = Petropoulos_metadata)
Petropoulos_TPM[['percent.mt']]<-PercentageFeatureSet(Petropoulos_TPM, pattern="^MT")

VlnPlot(Petropoulos_TPM, features=c('nFeature_RNA_tpm', 'nCount_RNA_tpm', 'percent.mt'))

DefaultAssay(Petropoulos)<-'RNA_qumi'

Petropoulos[['percent.mt']]<-PercentageFeatureSet(Petropoulos, pattern="^MT")

VlnPlot(Petropoulos, features=c('nFeature_RNA_qumi', 'nCount_RNA_qumi', 'percent.mt'))

Petropoulos<-NormalizeData(Petropoulos, normalization.method = "RC", scale.factor = 10e6)
Petropoulos<-FindVariableFeatures(Petropoulos)
all.genes<-rownames(Petropoulos)
Petropoulos<-ScaleData(Petropoulos, features=all.genes)

Petropoulos<-RunPCA(Petropoulos, features=VariableFeatures(object=Petropoulos))
DimHeatmap(Petropoulos, dims=1:15, cells=500, balanced=TRUE)

ElbowPlot(Petropoulos, ndims=50)

Petropoulos<-FindNeighbors(Petropoulos, dims=1:20)
Petropoulos<-FindClusters(Petropoulos, resolution=0.1)

Petropoulos<-RunUMAP(Petropoulos, dims=1:20)
DimPlot(Petropoulos, reduction="umap")
DimPlot(Petropoulos, group.by = "Age")

FeaturePlot(Petropoulos, features=c("SOX2", "SOX17", "GATA3"))

FeaturePlot(Petropoulos, features=c("NANOG", "SOX2", "POU5F1","TBXT",
                              "GATA4", "SOX17", "GATA6", "PDGFRA",
                              "GATA3", "GATA2", "TFAP2C", "CDX2"),
            order=TRUE, slot='data', min.cutoff = 'q25')

saveRDS(Petropoulos, "Petropoulos_qUMI.RDS")
