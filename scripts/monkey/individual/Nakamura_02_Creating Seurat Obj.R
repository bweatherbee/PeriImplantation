library(tidyverse)
library(Seurat)
library(readr)
qumi_matrix <- readRDS("C:/Users/baile/Desktop/realigning_for_int/CynMonkey/Nakamura/qumi_matrix.rds")
tpm_matrix <- readRDS("C:/Users/baile/Desktop/realigning_for_int/CynMonkey/Nakamura/tpm_matrix.rds")
metadata <- read_csv("meta.csv")

metadata<-as.data.frame(metadata[1:390,])
rownames(metadata)<-metadata$SampleID

Nakamura<-CreateSeuratObject(counts=qumi_matrix, project='Nakamura', assay='RNA',
                             meta.data=metadata)

Nakamura$percent.mt<-PercentageFeatureSet(Nakamura, pattern="^MT")
VlnPlot(Nakamura, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))

Nakamura<-NormalizeData(Nakamura, normalization.method = "RC", scale.factor = 10e6)
Nakamura<-FindVariableFeatures(Nakamura)
Nakamura<-ScaleData(Nakamura, features=rownames(Nakamura))
Nakamura<-RunPCA(Nakamura, features=VariableFeatures(Nakamura))
ElbowPlot(Nakamura)

Nakamura<-FindNeighbors(Nakamura, dims=1:10)
Nakamura<-FindClusters(Nakamura, resolution=0.2)

Nakamura<-RunUMAP(Nakamura, dims=1:5)
DimPlot(Nakamura, group.by = "Age")
DimPlot(Nakamura, group.by="cell_type")
DimPlot(Nakamura)
FeaturePlot(Nakamura, slot='data', features=c("SOX2", "NANOG", "POU5F1",
                                              "GATA3", "TFAP2C", "TFAP2A", "ISL1",
                                              "T", "MIXL1", "FOXH1", "EOMES",
                                              "GATA6", "GATA4", "SOX17", "COL6A1"), order=TRUE,max.cutoff = 'q98')

Idents(Nakamura)<-Nakamura$cell_type

saveRDS(Nakamura, 'Nakamura.RDS')
