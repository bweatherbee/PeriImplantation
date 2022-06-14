library(tidyverse)
library(Matrix)
library(Seurat)

t2g<-read_tsv("Cyn5.0_t2g.txt", col_names = c("RNAID", "GeneID", "Gene_name"))
genenamemap<-t2g[,2:3]
genenamemap<-unique(genenamemap)

MF111<-ReadMtx(mtx="./MFI111/counts_unfiltered/cells_x_genes.mtx", features="./MFI111/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI111/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF111<-t(as.matrix(MF111))
rownames(MF111)<-genenamemap$Gene_name
MF111<-CreateSeuratObject(counts=MF111, project="Ma", assay="RNA")
MF111$embryo<-"MFI11.1"
MF111$paper<-"Ma"
MF111$species<-"CynMonkey"
MF111$Age<-"D11"
MF111$stage<-"early post"
MF111$percent.mt<-PercentageFeatureSet(MF111, pattern="^MT")
VlnPlot(MF111, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF111<-subset(MF111, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF111, "MFI111.RDS")


MF112<-ReadMtx(mtx="./MFI112/counts_unfiltered/cells_x_genes.mtx", features="./MFI112/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI112/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF112<-t(as.matrix(MF112))
rownames(MF112)<-genenamemap$Gene_name
MF112<-CreateSeuratObject(counts=MF112, project="Ma", assay="RNA")
MF112$embryo<-"MFI11.2"
MF112$paper<-"Ma"
MF112$species<-"CynMonkey"
MF112$Age<-"D11"
MF112$stage<-"early post"
MF112$percent.mt<-PercentageFeatureSet(MF112, pattern="^MT")
VlnPlot(MF112, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF112<-subset(MF112, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF112, "MFI112.RDS")


MF121<-ReadMtx(mtx="./MFI121/counts_unfiltered/cells_x_genes.mtx", features="./MFI121/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI121/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF121<-t(as.matrix(MF121))
rownames(MF121)<-genenamemap$Gene_name
MF121<-CreateSeuratObject(counts=MF121, project="Ma", assay="RNA")
MF121$embryo<-"MFI12.1"
MF121$paper<-"Ma"
MF121$species<-"CynMonkey"
MF121$Age<-"D12"
MF121$stage<-"early post"
MF121$percent.mt<-PercentageFeatureSet(MF121, pattern="^MT")
VlnPlot(MF121, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF121<-subset(MF121, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF121, "MFI121.RDS")


MF122<-ReadMtx(mtx="./MFI122/counts_unfiltered/cells_x_genes.mtx", features="./MFI122/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI122/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF122<-t(as.matrix(MF122))
rownames(MF122)<-genenamemap$Gene_name
MF122<-CreateSeuratObject(counts=MF122, project="Ma", assay="RNA")
MF122$embryo<-"MFI12.2"
MF122$paper<-"Ma"
MF122$species<-"CynMonkey"
MF122$Age<-"D12"
MF122$stage<-"early post"
MF122$percent.mt<-PercentageFeatureSet(MF122, pattern="^MT")
VlnPlot(MF122, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF122<-subset(MF122, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF122, "MFI122.RDS")


MF123<-ReadMtx(mtx="./MFI123/counts_unfiltered/cells_x_genes.mtx", features="./MFI123/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI123/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF123<-t(as.matrix(MF123))
rownames(MF123)<-genenamemap$Gene_name
MF123<-CreateSeuratObject(counts=MF123, project="Ma", assay="RNA")
MF123$embryo<-"MFI12.3"
MF123$paper<-"Ma"
MF123$species<-"CynMonkey"
MF123$Age<-"D12"
MF123$stage<-"early post"
MF123$percent.mt<-PercentageFeatureSet(MF123, pattern="^MT")
VlnPlot(MF123, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF123<-subset(MF123, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF123, "MFI123.RDS")


MF124<-ReadMtx(mtx="./MFI124/counts_unfiltered/cells_x_genes.mtx", features="./MFI124/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI124/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF124<-t(as.matrix(MF124))
rownames(MF124)<-genenamemap$Gene_name
MF124<-CreateSeuratObject(counts=MF124, project="Ma", assay="RNA")
MF124$embryo<-"MFI12.4"
MF124$paper<-"Ma"
MF124$species<-"CynMonkey"
MF124$Age<-"D12"
MF124$stage<-"early post"
MF124$percent.mt<-PercentageFeatureSet(MF124, pattern="^MT")
VlnPlot(MF124, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF124<-subset(MF124, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF124, "MFI124.RDS")


MF125<-ReadMtx(mtx="./MFI125/counts_unfiltered/cells_x_genes.mtx", features="./MFI125/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI125/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF125<-t(as.matrix(MF125))
rownames(MF125)<-genenamemap$Gene_name
MF125<-CreateSeuratObject(counts=MF125, project="Ma", assay="RNA")
MF125$embryo<-"MFI12.5"
MF125$paper<-"Ma"
MF125$species<-"CynMonkey"
MF125$Age<-"D12"
MF125$stage<-"early post"
MF125$percent.mt<-PercentageFeatureSet(MF125, pattern="^MT")
VlnPlot(MF125, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF125<-subset(MF125, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF125, "MFI125.RDS")


MF131<-ReadMtx(mtx="./MFI131/counts_unfiltered/cells_x_genes.mtx", features="./MFI131/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI131/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF131<-t(as.matrix(MF131))
rownames(MF131)<-genenamemap$Gene_name
MF131<-CreateSeuratObject(counts=MF131, project="Ma", assay="RNA")
MF131$embryo<-"MFI13.1"
MF131$paper<-"Ma"
MF131$species<-"CynMonkey"
MF131$Age<-"D13"
MF131$stage<-"late post"
MF131$percent.mt<-PercentageFeatureSet(MF131, pattern="^MT")
VlnPlot(MF131, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF131<-subset(MF131, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF131, "MFI131.RDS")


MF132<-ReadMtx(mtx="./MFI132/counts_unfiltered/cells_x_genes.mtx", features="./MFI132/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI132/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF132<-t(as.matrix(MF132))
rownames(MF132)<-genenamemap$Gene_name
MF132<-CreateSeuratObject(counts=MF132, project="Ma", assay="RNA")
MF132$embryo<-"MFI13.2"
MF132$paper<-"Ma"
MF132$species<-"CynMonkey"
MF132$Age<-"D13"
MF132$stage<-"late post"
MF132$percent.mt<-PercentageFeatureSet(MF132, pattern="^MT")
VlnPlot(MF132, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF132<-subset(MF132, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF132, "MFI132.RDS")


MF133<-ReadMtx(mtx="./MFI133/counts_unfiltered/cells_x_genes.mtx", features="./MFI133/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI133/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF133<-t(as.matrix(MF133))
rownames(MF133)<-genenamemap$Gene_name
MF133<-CreateSeuratObject(counts=MF133, project="Ma", assay="RNA")
MF133$embryo<-"MFI13.3"
MF133$paper<-"Ma"
MF133$species<-"CynMonkey"
MF133$Age<-"D13"
MF133$stage<-"late post"
MF133$percent.mt<-PercentageFeatureSet(MF133, pattern="^MT")
VlnPlot(MF133, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF133<-subset(MF133, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF133, "MFI133.RDS")


MF141<-ReadMtx(mtx="./MFI141/counts_unfiltered/cells_x_genes.mtx", features="./MFI141/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI141/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF141<-t(as.matrix(MF141))
rownames(MF141)<-genenamemap$Gene_name
MF141<-CreateSeuratObject(counts=MF141, project="Ma", assay="RNA")
MF141$embryo<-"MFI14.1"
MF141$paper<-"Ma"
MF141$species<-"CynMonkey"
MF141$Age<-"D14"
MF141$stage<-"late post"
MF141$percent.mt<-PercentageFeatureSet(MF141, pattern="^MT")
VlnPlot(MF141, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF141<-subset(MF141, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF141, "MFI141.RDS")


MF142<-ReadMtx(mtx="./MFI142/counts_unfiltered/cells_x_genes.mtx", features="./MFI142/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI142/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF142<-t(as.matrix(MF142))
rownames(MF142)<-genenamemap$Gene_name
MF142<-CreateSeuratObject(counts=MF142, project="Ma", assay="RNA")
MF142$embryo<-"MFI14.2"
MF142$paper<-"Ma"
MF142$species<-"CynMonkey"
MF142$Age<-"D14"
MF142$stage<-"late post"
MF142$percent.mt<-PercentageFeatureSet(MF142, pattern="^MT")
VlnPlot(MF142, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF142<-subset(MF142, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF142, "MFI142.RDS")


MF143<-ReadMtx(mtx="./MFI143/counts_unfiltered/cells_x_genes.mtx", features="./MFI143/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI143/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF143<-t(as.matrix(MF143))
rownames(MF143)<-genenamemap$Gene_name
MF143<-CreateSeuratObject(counts=MF143, project="Ma", assay="RNA")
MF143$embryo<-"MFI14.3"
MF143$paper<-"Ma"
MF143$species<-"CynMonkey"
MF143$Age<-"D14"
MF143$stage<-"late post"
MF143$percent.mt<-PercentageFeatureSet(MF143, pattern="^MT")
VlnPlot(MF143, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF143<-subset(MF143, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF143, "MFI143.RDS")


MF144<-ReadMtx(mtx="./MFI144/counts_unfiltered/cells_x_genes.mtx", features="./MFI144/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI144/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF144<-t(as.matrix(MF144))
rownames(MF144)<-genenamemap$Gene_name
MF144<-CreateSeuratObject(counts=MF144, project="Ma", assay="RNA")
MF144$embryo<-"MFI14.4"
MF144$paper<-"Ma"
MF144$species<-"CynMonkey"
MF144$Age<-"D14"
MF144$stage<-"late post"
MF144$percent.mt<-PercentageFeatureSet(MF144, pattern="^MT")
VlnPlot(MF144, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF144<-subset(MF144, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF144, "MFI144.RDS")


MF161<-ReadMtx(mtx="./MFI161/counts_unfiltered/cells_x_genes.mtx", features="./MFI161/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI161/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF161<-t(as.matrix(MF161))
rownames(MF161)<-genenamemap$Gene_name
MF161<-CreateSeuratObject(counts=MF161, project="Ma", assay="RNA")
MF161$embryo<-"MFI16.1"
MF161$paper<-"Ma"
MF161$species<-"CynMonkey"
MF161$Age<-"D16"
MF161$stage<-"gast"
MF161$percent.mt<-PercentageFeatureSet(MF161, pattern="^MT")
VlnPlot(MF161, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF161<-subset(MF161, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF161, "MFI161.RDS")


MF162<-ReadMtx(mtx="./MFI162/counts_unfiltered/cells_x_genes.mtx", features="./MFI162/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI162/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF162<-t(as.matrix(MF162))
rownames(MF162)<-genenamemap$Gene_name
MF162<-CreateSeuratObject(counts=MF162, project="Ma", assay="RNA")
MF162$embryo<-"MFI16.2"
MF162$paper<-"Ma"
MF162$species<-"CynMonkey"
MF162$Age<-"D16"
MF162$stage<-"gast"
MF162$percent.mt<-PercentageFeatureSet(MF162, pattern="^MT")
VlnPlot(MF162, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF162<-subset(MF162, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF162, "MFI162.RDS")


MF163<-ReadMtx(mtx="./MFI163/counts_unfiltered/cells_x_genes.mtx", features="./MFI163/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI163/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF163<-t(as.matrix(MF163))
rownames(MF163)<-genenamemap$Gene_name
MF163<-CreateSeuratObject(counts=MF163, project="Ma", assay="RNA")
MF163$embryo<-"MFI16.3"
MF163$paper<-"Ma"
MF163$species<-"CynMonkey"
MF163$Age<-"D16"
MF163$stage<-"gast"
MF163$percent.mt<-PercentageFeatureSet(MF163, pattern="^MT")
VlnPlot(MF163, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF163<-subset(MF163, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF163, "MFI163.RDS")


MF164<-ReadMtx(mtx="./MFI164/counts_unfiltered/cells_x_genes.mtx", features="./MFI164/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI164/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF164<-t(as.matrix(MF164))
rownames(MF164)<-genenamemap$Gene_name
MF164<-CreateSeuratObject(counts=MF164, project="Ma", assay="RNA")
MF164$embryo<-"MFI164"
MF164$paper<-"Ma"
MF164$species<-"CynMonkey"
MF164$Age<-"D16"
MF164$stage<-"gast"
MF164$percent.mt<-PercentageFeatureSet(MF164, pattern="^MT")
VlnPlot(MF164, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF164<-subset(MF164, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF164, "MFI164.RDS")


MF171<-ReadMtx(mtx="./MFI171/counts_unfiltered/cells_x_genes.mtx", features="./MFI171/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI171/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF171<-t(as.matrix(MF171))
rownames(MF171)<-genenamemap$Gene_name
MF171<-CreateSeuratObject(counts=MF171, project="Ma", assay="RNA")
MF171$embryo<-"MFI17.1"
MF171$paper<-"Ma"
MF171$species<-"CynMonkey"
MF171$Age<-"D17"
MF171$stage<-"gast"
MF171$percent.mt<-PercentageFeatureSet(MF171, pattern="^MT")
VlnPlot(MF171, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF171<-subset(MF171, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF171, "MFI171.RDS")


MF172<-ReadMtx(mtx="./MFI172/counts_unfiltered/cells_x_genes.mtx", features="./MFI172/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI172/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF172<-t(as.matrix(MF172))
rownames(MF172)<-genenamemap$Gene_name
MF172<-CreateSeuratObject(counts=MF172, project="Ma", assay="RNA")
MF172$embryo<-"MFI17.2"
MF172$paper<-"Ma"
MF172$species<-"CynMonkey"
MF172$Age<-"D17"
MF172$stage<-"gast"
MF172$percent.mt<-PercentageFeatureSet(MF172, pattern="^MT")
VlnPlot(MF172, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF172<-subset(MF172, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF172, "MFI172.RDS")


MF173<-ReadMtx(mtx="./MFI173/counts_unfiltered/cells_x_genes.mtx", features="./MFI173/counts_unfiltered/cells_x_genes.barcodes.txt",
               cells="./MFI173/counts_unfiltered/cells_x_genes.genes.txt",
               feature.column=1)
MF173<-t(as.matrix(MF173))
rownames(MF173)<-genenamemap$Gene_name
MF173<-CreateSeuratObject(counts=MF173, project="Ma", assay="RNA")
MF173$embryo<-"MFI17.3"
MF173$paper<-"Ma"
MF173$species<-"CynMonkey"
MF173$Age<-"D17"
MF173$stage<-"gast"
MF173$percent.mt<-PercentageFeatureSet(MF173, pattern="^MT")
VlnPlot(MF173, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
MF173<-subset(MF173, subset=nFeature_RNA>1000 & percent.mt<20)
saveRDS(MF173, "MFI173.RDS")

Ma<-merge(x=MF111, y=c(MF112, MF121, MF122, MF123, MF124, MF125, MF131,
                       MF132, MF133, MF141, MF142, MF143, MF144, MF161,
                       MF162, MF163, MF164, MF171, MF172, MF173))

rm(MF111, MF112, MF121, MF122, MF123, MF124, MF125, 
   MF131, MF132, MF133, MF141, MF142, MF143, MF144, 
   MF161, MF162, MF163, MF164, MF171, MF172, MF173)

Ma<-subset(Ma, subset=nFeature_RNA>2000 & nCount_RNA<25000)

VlnPlot(Ma, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Ma<-NormalizeData(Ma, normalization.method = "RC", scale.factor = 10e6)
Ma<-FindVariableFeatures(Ma)
Ma<-ScaleData(Ma, features=rownames(Ma))
Ma<-RunPCA(Ma, features=VariableFeatures(Ma))
ElbowPlot(Ma)

Ma<-FindNeighbors(Ma, dims=1:10)
Ma<-FindClusters(Ma, resolution=0.1)

Ma<-RunUMAP(Ma, dims=1:20)
DimPlot(Ma, group.by = "Age")
DimPlot(Ma, group.by="embryo")
DimPlot(Ma)
FeaturePlot(Ma, slot='data', features=c("SOX2", "NANOG", "POU5F1",
                                          "GATA3", "TFAP2C", "TFAP2A", "ISL1",
                                          "T", "MIXL1", "FOXH1", "EOMES",
                                          "GATA6", "GATA4", "SOX17", "COL6A1"), order=TRUE,max.cutoff = 'q98')

new.cluster.ids<-c("EPI", "EXMC", "TE", "HYPO")
names(new.cluster.ids)<-levels(Ma)
Ma<-RenameIdents(Ma, new.cluster.ids)
DimPlot(Ma)
Ma$cell_type<-Idents(Ma)
saveRDS(Ma, "Ma.RDS")
