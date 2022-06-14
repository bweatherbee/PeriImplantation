library(tidyverse)
library(Matrix)
library(Seurat)

setwd("C:/Users/baile/Desktop/realigning_for_int/CynMonkey")
t2g<-read_tsv("Cyn5.0_t2g.txt", col_names = c("RNAID", "GeneID", "Gene_name"))
genenamemap<-t2g[,2:3]
genenamemap<-unique(genenamemap)

D10_B1<-ReadMtx(mtx="./Yang/WT_d10_B1/counts_unfiltered/cells_x_genes.mtx", features="./Yang/WT_d10_B1/counts_unfiltered/cells_x_genes.barcodes.txt",
                cells="./Yang/WT_d10_B1/counts_unfiltered/cells_x_genes.genes.txt",
                feature.column=1)
D10_B1<-t((D10_B1))
D10_B1@Dimnames[[1]]<-genenamemap$Gene_name
D10_B1<-CreateSeuratObject(counts=D10_B1, project="Yang", assay = "RNA")
D10_B1$embryo<-"D10_B1"
D10_B1$paper<-"Yang"
D10_B1$species<-"CynMonkey"
D10_B1$Age<-"D10"
D10_B1$stage<-"early post"
D10_B1$percent.mt<-PercentageFeatureSet(D10_B1, pattern="^MT")
max(D10_B1$percent.mt, na.rm=T)
max(D10_B1$nFeature_RNA)
D10_B1<-subset(D10_B1, subset=nFeature_RNA>1000 & percent.mt<20)
VlnPlot(D10_B1, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
saveRDS(D10_B1, "./Yang/D10_B1.rds")



D10_B2<-ReadMtx(mtx="./Yang/WT_D10_B2/counts_unfiltered/cells_x_genes.mtx", features="./Yang/WT_D10_B2/counts_unfiltered/cells_x_genes.barcodes.txt",
                cells="./Yang/WT_D10_B2/counts_unfiltered/cells_x_genes.genes.txt",
                feature.column=1)
D10_B2<-t((D10_B2))
D10_B2@Dimnames[[1]]<-genenamemap$Gene_name
D10_B2<-CreateSeuratObject(counts=D10_B2, project="Yang", assay = "RNA")
D10_B2$embryo<-"D10_B2"
D10_B2$paper<-"Yang"
D10_B2$species<-"CynMonkey"
D10_B2$Age<-"D10"
D10_B2$stage<-"early post"
D10_B2$percent.mt<-PercentageFeatureSet(D10_B2, pattern="^MT")
max(D10_B2$percent.mt, na.rm=T)
max(D10_B2$nFeature_RNA)
D10_B2<-subset(D10_B2, subset=nFeature_RNA>1000 & percent.mt<20)
VlnPlot(D10_B2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
saveRDS(D10_B2, "./Yang/D10_B2.rds")

D12_B1<-ReadMtx(mtx="./Yang/WT_D12_B1/counts_unfiltered/cells_x_genes.mtx", features="./Yang/WT_D12_B1/counts_unfiltered/cells_x_genes.barcodes.txt",
                cells="./Yang/WT_D12_B1/counts_unfiltered/cells_x_genes.genes.txt",
                feature.column=1)
D12_B1<-t((D12_B1))
D12_B1@Dimnames[[1]]<-genenamemap$Gene_name
D12_B1<-CreateSeuratObject(counts=D12_B1, project="Yang", assay = "RNA")
D12_B1$embryo<-"D12_B1"
D12_B1$paper<-"Yang"
D12_B1$species<-"CynMonkey"
D12_B1$Age<-"D12"
D12_B1$stage<-"early post"
D12_B1$percent.mt<-PercentageFeatureSet(D12_B1, pattern="^MT")
max(D12_B1$percent.mt, na.rm=T)
max(D12_B1$nFeature_RNA)
D12_B1<-subset(D12_B1, subset=nFeature_RNA>1000 & percent.mt<20)
VlnPlot(D12_B1, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
saveRDS(D12_B1, "./Yang/D12_B1.rds")


D12_B2<-ReadMtx(mtx="./Yang/WT_D12_B2/counts_unfiltered/cells_x_genes.mtx", features="./Yang/WT_D12_B2/counts_unfiltered/cells_x_genes.barcodes.txt",
                cells="./Yang/WT_D12_B2/counts_unfiltered/cells_x_genes.genes.txt",
                feature.column=1)
D12_B2<-t((D12_B2))
D12_B2@Dimnames[[1]]<-genenamemap$Gene_name
D12_B2<-CreateSeuratObject(counts=D12_B2, project="Yang", assay = "RNA")
D12_B2$embryo<-"D12_B2"
D12_B2$paper<-"Yang"
D12_B2$species<-"CynMonkey"
D12_B2$Age<-"D12"
D12_B2$stage<-"early post"
D12_B2$percent.mt<-PercentageFeatureSet(D12_B2, pattern="^MT")
max(D12_B2$percent.mt, na.rm=T)
max(D12_B2$nFeature_RNA)
D12_B2<-subset(D12_B2, subset=nFeature_RNA>1000 & percent.mt<20)
VlnPlot(D12_B2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
saveRDS(D12_B2, "./Yang/D12_B2.rds")


D14_B1<-ReadMtx(mtx="./Yang/WT_D14_B1/counts_unfiltered/cells_x_genes.mtx", features="./Yang/WT_D14_B1/counts_unfiltered/cells_x_genes.barcodes.txt",
                cells="./Yang/WT_D14_B1/counts_unfiltered/cells_x_genes.genes.txt",
                feature.column=1)
D14_B1<-t((D14_B1))
D14_B1@Dimnames[[1]]<-genenamemap$Gene_name
D14_B1<-CreateSeuratObject(counts=D14_B1, project="Yang", assay = "RNA")
D14_B1$embryo<-"D14_B1"
D14_B1$paper<-"Yang"
D14_B1$species<-"CynMonkey"
D14_B1$Age<-"D14"
D14_B1$stage<-"late post"
D14_B1$percent.mt<-PercentageFeatureSet(D14_B1, pattern="^MT")
max(D14_B1$percent.mt, na.rm=T)
max(D14_B1$nFeature_RNA)
D14_B1<-subset(D14_B1, subset=nFeature_RNA>1000 & percent.mt<20)
VlnPlot(D14_B1, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
saveRDS(D14_B1, "./Yang/D14_B1.rds")


D14_B2<-ReadMtx(mtx="./Yang/WT_D14_B2/counts_unfiltered/cells_x_genes.mtx", features="./Yang/WT_D14_B2/counts_unfiltered/cells_x_genes.barcodes.txt",
                cells="./Yang/WT_D14_B2/counts_unfiltered/cells_x_genes.genes.txt",
                feature.column=1)
D14_B2<-t((D14_B2))
D14_B2@Dimnames[[1]]<-genenamemap$Gene_name
D14_B2<-CreateSeuratObject(counts=D14_B2, project="Yang", assay = "RNA")
D14_B2$embryo<-"D14_B2"
D14_B2$paper<-"Yang"
D14_B2$species<-"CynMonkey"
D14_B2$Age<-"D14"
D14_B2$stage<-"late post"
D14_B2$percent.mt<-PercentageFeatureSet(D14_B2, pattern="^MT")
max(D14_B2$percent.mt, na.rm=T)
max(D14_B2$nFeature_RNA)
D14_B2<-subset(D14_B2, subset=nFeature_RNA>1000 & percent.mt<20)
VlnPlot(D14_B2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
saveRDS(D14_B2, "./Yang/D14_B2.rds")


Yang<-merge(D10_B1, y=c(D10_B2, D12_B1, D12_B2, D14_B1, D14_B2))

rm(D10_B1, D10_B2, D12_B1, D12_B2, D14_B1, D14_B2)

VlnPlot(Yang, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Yang<-NormalizeData(Yang, normalization.method = "RC", scale.factor = 10e6)
Yang<-FindVariableFeatures(Yang)
Yang<-ScaleData(Yang, features=rownames(Yang))
Yang<-RunPCA(Yang, features=VariableFeatures(Yang))
ElbowPlot(Yang)

Yang<-FindNeighbors(Yang, dims=1:15)
Yang<-FindClusters(Yang, resolution=0.02)

Yang<-RunUMAP(Yang, dims=1:20)
DimPlot(Yang, group.by = "Age")
DimPlot(Yang, group.by="embryo")
DimPlot(Yang)
FeaturePlot(Yang, slot='data', features=c("SOX2", "NANOG", "POU5F1",
                                        "GATA3", "TFAP2C", "TFAP2A", "ISL1",
                                        "T", "MIXL1", "FOXH1", "EOMES",
                                        "GATA6", "GATA4", "SOX17", "COL6A1"), order=TRUE,max.cutoff = 'q98')

new.cluster.ids<-c("TE", "EPI", "EXMC", "HYPO")
names(new.cluster.ids)<-levels(Yang)
Yang<-RenameIdents(Yang, new.cluster.ids)
DimPlot(Yang)
Yang$cell_type<-Idents(Yang)
saveRDS(Yang, "./Yang/Yang.RDS")
