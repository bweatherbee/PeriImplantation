library(Seurat)
library(biomaRt)
library(org.Hs.eg.db)
#BiocManager::install("rWikiPathways")
library(rWikiPathways)

merged <- readRDS("D:/realigning_for_int/FOR 2022 SIGNALING PAPER/human/merged_integrated.RDS")
hs.pathways <- listPathways('Homo sapiens')


ensembl=useMart("ensembl")
datasets<-listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl", mart=ensembl)

AMPK<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP1403', 'En'), mart=ensembl)
AMPK<-list(unique(AMPK$external_gene_name))

BMP1<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP1425', 'En'), mart=ensembl)
BMP1<-list(unique(BMP1$external_gene_name))
BMP2<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP2760', 'En'), mart=ensembl)
BMP2<-list(unique(BMP2$external_gene_name))

NODAL<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP2675', 'En'), mart=ensembl)
NODAL<-list(unique(NODAL$external_gene_name))
Activin<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP2791', 'En'), mart=ensembl)
Activin<-list(unique(Activin$external_gene_name))

TGFBR_C<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP2742', 'En'), mart=ensembl)
TGFBR_C<-list(unique(TGFBR_C$external_gene_name))
TGFB<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP366', 'En'), mart=ensembl)
TGFB<-list(unique(TGFB$external_gene_name))
TGFBR<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                 filters="ensembl_gene_id",
                 values=getXrefList('WP560', 'En'), mart=ensembl)
TGFBR<-list(unique(TGFBR$external_gene_name))

IGF1R<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                   filters="ensembl_gene_id",
                   values=getXrefList('WP2677', 'En'), mart=ensembl)
IGF1R<-list(unique(IGF1R$external_gene_name))

Notch<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                   filters="ensembl_gene_id",
                   values=getXrefList('WP268', 'En'), mart=ensembl)
Notch<-list(unique(Notch$external_gene_name))

Hippo<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                   filters="ensembl_gene_id",
                   values=getXrefList('WP2714', 'En'), mart=ensembl)
Hippo<-list(unique(Hippo$external_gene_name))
Hippo_Yap<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                       filters="ensembl_gene_id",
                       values=getXrefList('WP4537', 'En'), mart=ensembl)
Hippo_Yap<-list(unique(Hippo_Yap$external_gene_name))

Hippo_reg<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP4540', 'En'), mart=ensembl)
Hippo_reg<-list(unique(Hippo_reg$external_gene_name))

MAPK<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                  filters="ensembl_gene_id",
                  values=getXrefList('WP382', 'En'), mart=ensembl)
MAPK<-list(unique(MAPK$external_gene_name))
MAPK4_6<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                  filters="ensembl_gene_id",
                  values=getXrefList('WP3307', 'En'), mart=ensembl)
MAPK4_6<-list(unique(MAPK4_6$external_gene_name))

FGFR1<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP3335', 'En'), mart=ensembl)
FGFR1<-list(unique(FGFR1$external_gene_name))
FGFR2<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
             filters="ensembl_gene_id",
             values=getXrefList('WP3338', 'En'), mart=ensembl)
FGFR2<-list(unique(FGFR2$external_gene_name))
FGFR3<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
             filters="ensembl_gene_id",
             values=getXrefList('WP3336', 'En'), mart=ensembl)
FGFR3<-list(unique(FGFR3$external_gene_name))
FGFR4<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
             filters="ensembl_gene_id",
             values=getXrefList('WP3334', 'En'), mart=ensembl)
FGFR4<-list(unique(FGFR4$external_gene_name))

PI3K_AKT<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                     filters="ensembl_gene_id",
                     values=getXrefList('WP4172', 'En'), mart=ensembl)
PI3K_AKT<-list(unique(PI3K_AKT$external_gene_name))

mTOR<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                 filters="ensembl_gene_id",
                 values=getXrefList('WP3318', 'En'), mart=ensembl)
mTOR<-list(unique(mTOR$external_gene_name))

Hh<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                 filters="ensembl_gene_id",
                 values=getXrefList('WP4249', 'En'), mart=ensembl)
Hh<-list(unique(Hh$external_gene_name))

RA<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
             filters="ensembl_gene_id",
             values=getXrefList('WP4421', 'En'), mart=ensembl)
RA<-list(unique(RA$external_gene_name))

WNT1<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
             filters="ensembl_gene_id",
             values=getXrefList('WP363', 'En'), mart=ensembl)
WNT1<-list(unique(WNT1$external_gene_name))
WNT2<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP428', 'En'), mart=ensembl)
WNT2<-list(unique(WNT2$external_gene_name))
WNT_bcat_ind<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP3558', 'En'), mart=ensembl)
WNT_bcat_ind<-list(unique(WNT_bcat_ind$external_gene_name))
WNT_TCF_dep<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="ensembl_gene_id",
            values=getXrefList('WP3339', 'En'), mart=ensembl)
WNT_TCF_dep<-list(unique(WNT_TCF_dep$external_gene_name))



Pathways<-list("AMPK"=AMPK[[1]], "BMP1"=BMP1[[1]], "BMP2"=BMP2[[1]], "NODAL"=NODAL[[1]], 
               "Activin"=Activin[[1]], "TGFBR_C"=TGFBR_C[[1]], "TGFB"=TGFB[[1]], 
               "TGFBR"=TGFBR[[1]], "IGF1R"=IGF1R[[1]], "Notch"=Notch[[1]], "Hippo"=Hippo[[1]], 
               "Hippo_Yap"=Hippo_Yap[[1]], "Hippo_reg"=Hippo_reg[[1]], "MAPK"= MAPK[[1]], 
               "MAPK4_6"=MAPK4_6[[1]], "FGFR1"=FGFR1[[1]], "FGFR2"=FGFR2[[1]], "FGFR3"=FGFR3[[1]], 
               "FGFR4"=FGFR4[[1]], "PI3K_AKT"=PI3K_AKT[[1]],"mTOR"=mTOR[[1]], "Hh"=Hh[[1]], 
               "RA"=RA[[1]], "WNT1"=WNT1[[1]], "WNT2"=WNT2[[1]], "WNT_bcat_ind"=WNT_bcat_ind[[1]], 
               "WNT_TCF_dep"=WNT_TCF_dep[[1]])

DefaultAssay(merged)<-"SCT"
merged_scored<-AddModuleScore(merged, features=Pathways, assay="SCT", name="signaling")
#options(scipen = 999)
signaling_matrix<-merged_scored@meta.data[,24:50]
colnames(signaling_matrix)<-c("AMPK", "BMP1", "BMP2", "NODAL", "Activin", "TGFBR_C", "TGFB", "TGFBR",
                              "IGF1R", "Notch", "Hippo", "Hippo_Yap", "Hippo_reg", "MAPK", "MAPK4_6",
                              "FGFR1", "FGFR2", "FGFR3", "FGFR4", "PI3K_AKT", "mTOR", "Hh", "RA",
                              "WNT1", "WNT2", "WNT_bcat_ind", "WNT_TCF_dep")
merged_scored<-AddMetaData(merged, metadata=signaling_matrix)
saveRDS(signaling_matrix, 'wikipathways_modulescores.rds')


library(Scillus)
library(RColorBrewer)
merged_scored@assays[['RNA']]<-merged_scored@assays$UMI_qumi
plot_measure(dataset = merged_scored, 
             measures = c("AMPK", "BMP1", "BMP2", "NODAL", "Activin", "TGFBR_C", "TGFB", "TGFBR",
                          "IGF1R", "Notch", "Hippo", "Hippo_Yap", "Hippo_reg", "MAPK", "MAPK4_6",
                          "FGFR1", "FGFR2", "FGFR3", "FGFR4", "PI3K_AKT", "mTOR", "Hh", "RA",
                          "WNT1", "WNT2", "WNT_bcat_ind", "WNT_TCF_dep"), 
             group_by = "new_cell_type",
             pal_setup = "Dark2")

plot_measure(dataset = merged_scored, 
             measures = c("AMPK"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("BMP1"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1") +theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("BMP2"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("NODAL"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("Activin"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("TGFBR_C"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("TGFB"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("TGFBR"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("IGF1R"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("Notch"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("Hippo"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("Hippo_Yap"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("Hippo_reg"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("MAPK"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("MAPK4_6"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("FGFR1"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("FGFR2"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("FGFR3"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("FGFR4"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("PI3K_AKT"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("mTOR"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("Hh"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("RA"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("WNT1"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("WNT2"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("WNT_bcat_ind"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("WNT_TCF_dep"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")

FeaturePlot(merged_scored, features=c("AMPK", "BMP1", "BMP2", "NODAL", "Activin", "TGFBR_C", "TGFB", "TGFBR",
                               "IGF1R", "Notch", "Hippo", "Hippo_Yap", "Hippo_reg", "MAPK", "MAPK4_6",
                               "FGFR1", "FGFR2", "FGFR3", "FGFR4", "PI3K_AKT", "mTOR", "Hh", "RA",
                               "WNT1", "WNT2", "WNT_bcat_ind", "WNT_TCF_dep"),
            min.cutoff = 'q5', max.cutoff = 'q98', order=TRUE, cols=c("gray88", "green3"), pt.size=2)

Scored<-CreateSeuratObject(counts=t(signaling_matrix),assay='signaling', meta.data = merged@meta.data)

DefaultAssay(Scored)<-"signaling"
Scored<-ScaleData(Scored,assay='signaling')

plot_heatmap(dataset=Scored, markers=c("AMPK", "BMP1", "BMP2", "NODAL", "Activin", "TGFBR-C", "TGFB", "TGFBR",
                                              "IGF1R", "Notch", "Hippo", "Hippo-Yap", "Hippo-reg", "MAPK", "MAPK4-6",
                                              "FGFR1", "FGFR2", "FGFR3", "FGFR4", "PI3K-AKT", "mTOR", "Hh", "RA",
                                              "WNT1", "WNT2", "WNT-bcat-ind", "WNT-TCF-dep"),
             sort_var=c("new_cell_type", "stage", "paper"),
             anno_var=c("new_cell_type", "stage", "paper"),
             anno_colors = c("Set1", "Set2", "Set3"),
             hm_limit = c(-1,0,1))

Idents(Scored)<-Scored$new_cell_type
avg_activity<-as.data.frame(AverageExpression(Scored, assays="signaling", slot="counts", group.by = "new_cell_type", add.ident="Age"))
pheatmap::pheatmap(t(scale(t(avg_activity))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
avg_activity2<-as.data.frame(AverageExpression(Scored, assays="signaling", slot="counts", group.by="new_cell_type"))
pheatmap::pheatmap(t(scale(t(avg_activity2))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap::pheatmap(t(scale(t(avg_activity2[,1:3]))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

metadata<-merged@meta.data

Idents(Scored)<-Scored$new_cell_type

avg_activity3<-as.data.frame(AverageExpression(Scored, assays="signaling", slot="counts", group.by="new_cell_type", add.ident="stage"))

pheatmap::pheatmap(t(scale(t(avg_activity3))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap::pheatmap(t(scale(t(avg_activity3[,1:12]))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))



### Plotting groups of vlns for specific wikipathway module components ###
BMP1_genes<-c( "BMPR1A","RUNX2","BMP2", "SMAD6","BMPR1B","TOB1", "SMAD4","SMAD1","NOG","TOB2", "SMURF1","BMPR2")
BMP2_genes<-c("ZFYVE16","UBE2D1","SMAD7", "CHRDL1",
              "AMH","BMPR1A","SMURF2","UBE2D3","SMAD5",
              "ACVR2B","SMAD9", "ACVR2A","BMP2","AMHR2",
              "SMAD6", "BMPR1B","ACVRL1","SMAD4", "CER1","SKI",
              "BMP10","FSTL1", "SMAD1", "GREM2", "NOG","SMURF1","BMPR2", "GDF2")
NODAL_genes<-c("ACVR2B","FOXO3","ACVR2A","ACVR1C","GDF1", "ACVR1B","CFC1",
               "PCSK6","FURIN","SMAD4","LEFTY2","CER1", "NODAL","FOXH1",
               "SMAD3","SMAD2","DRAP1","DAND5","TDGF1","LEFTY1")
Activin_genes<-c("FSTL3", "ACVR2B","ACVR2A","INHBA","ACVR1C","FST","ACVR1B", "SMAD4", "FOXH1", "INHBB", "SMAD3", "SMAD2", "DRAP1")
TGFBR_C_genes<-c("STRAP","NEDD4L", "RHOA", "PRKCZ","XPO1", "PPP1R15A",
                 "FKBP1A", "BAMBI","SMAD7","PARD6A", "STUB1","ARHGEF18",
                 "TGFB1","TGFBR1", "MTMR4","SMURF2", "CBL","UCHL5", 
                 "PMEPA1", "NEDD8","UBE2M","USP15","FURIN","SMAD4", 
                 "CGN","RPS27A", "PARD3","UBC","F11R", "TGFBR2",
                 "SMAD3","UBB","PPP1CA", "SMAD2","PPP1CC", "SMURF1",
                 "PPP1CB","UBA52")
TGFB_genes<-c("CREBBP", "STRAP","PIAS1","MAP2K3", "ZFYVE16","TNC", 
              "NEDD4L", "FOXP3","MAPK9","BCAR1","CUL1", "MAP2K4",
              "KLF6", "RHOA", "ROCK1","MEF2A","TGFBR3", "CDC42", 
              "PIAS2","ITCH", "RBL1", "MEF2C","MAPK1","TAB1",
              "RBX1", "EP300","SNW1", "PPM1A","SMAD7","PARD6A",
              "AXIN1","UBE2I","RBL2", "MAP4K1", "TGFB1","PIK3R2",
              "CAV1", "MET","TGFBR1", "MAPK8","SMURF2", "MAP2K6",
              "CCND1","NEDD9","MAPK14", "SKP1", "PRKAR2A","SPTBN1",
              "FN1","SOS1", "ATF2", "SUMO1","HDAC1","UCHL5", 
              "COPS5","WWP1", "CDKN1A", "NUP153", "RUNX2","FOSB",
              "CITED1", "TRAP1","NUP214", "MAP2K2", "JUND", "RAF1",
              "ITGB4","E2F5", "ETS1", "MAP3K7", "RAC1", "SKIL",
              "MYC","YAP1", "THBS1","STAMBPL1","PML","TGFB1I1", 
              "TP53", "MAPK4","SMAD4","SIK1", "APP","AKT1",
              "PIK3R1", "CDKN2B", "ZEB1", "ITGB1","PDK1", "DAB2",
              "KLF10","ZFYVE9", "RNF111", "CCNB2","SKI","SHC1",
              "FOXH1","ATF3", "TGFBR2", "SNIP1","ITGA2","TERT",
              "COL1A2", "BTRC", "SMAD3","MAP2K1", "SIN3A","PTK2",
              "ZEB2", "CDK1", "FOS","JUNB", "KLF11","TRAF6", 
              "SMAD2","EID2", "TGIF1","JUN","GRB2", "PAK2",
              "PJA1", "LIMK2","HGS","SP1","MMP1", "SRC", 
              "TFDP1","SMURF1", "E2F4", "ITGB3","MMP12","DCP1A")
TGFBR_genes<-c("CREBBP", "RUNX3","LTBP1","MAPK9","TFE3", "TGFBR3",
               "FKBP1A", "BAMBI","EP300","SMAD7","MAPK3","ZNF423",
               "TGFB1","SERPINE1","TGFBR1", "ENG","NFKB1","IFNG",
               "SMAD5","ITGB6","STAT1","SPP1", "SMAD9","INHBA", 
               "RUNX2","WNT1", "BMP4", "LIF","FST","SKIL",
               "THBS1","SMAD6","LEF1", "EGF","SMAD4","LEFTY2",
               "ZFYVE9", "SKI","FOXH1","JAK1", "TGFBR2", "SMAD3", 
               "CTNNB1", "STAT3","ZEB2", "FOS","SMAD1","HRAS",
               "SMAD2","TGIF1","JUN","NOG","MIR302A","TNF", 
               "LEFTY1")
Notch_genes<-c("DVL2","CREBBP","NOTCH3","PSEN1","DLL3","DTX2","MFNG",
               "JAG1","NUMBL","LFNG","DVL1","KAT2A","DTX4","KAT2B",
               "HES1","HDAC1","APH1A","KCNJ5","RBPJL","DLL4","INPP5K",
               "NUMB","NOTCH2","DTX1","APH1B","PSEN2","NOTCH1","ADAM17",
               "CTBP1","MAML1","DVL3","NCSTN","DTX3L","RBPJ","RFNG",
               "PTCRA","CTBP2","DTX3","JAG2","NCOR2","HDAC2","MAML3",
               "HES5","DLL1","NOTCH4","TNF")



merged@assays[['RNA']]<-merged@assays$UMI_qumi
DefaultAssay(merged)<-'SCT'
plot_measure(dataset = merged, 
             measures = Activin_genes, 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")

plot_measure(dataset = merged, 
             measures = BMP1_genes, 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")

plot_measure(dataset = merged, 
             measures = BMP2_genes, 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")

plot_measure(dataset = merged, 
             measures = NODAL_genes, 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")

plot_measure(dataset = merged, 
             measures = Notch_genes, 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")

plot_measure(dataset = merged, 
             measures = TGFB_genes, 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")

plot_measure(dataset = merged, 
             measures = TGFBR_C_genes, 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")

plot_measure(dataset = merged, 
             measures = TGFBR_genes, 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")


### Plotting selected genes from above lists ###
merged@assays[['RNA']]<-merged@assays$UMI_qumi
DefaultAssay(merged)<-'RNA'
merged<-NormalizeData(merged)

#BMP-related

plot_measure(dataset = merged, 
             measures ="BMPR1A", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="BMPR1B", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="BMPR2", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="AMHR2", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="BMP2", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="BMP4", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="BMP6", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="SMAD1", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="SMAD5", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="SMAD9", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="ID2", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="NOG", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="GREM2", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()


#Broad TGFB Superfamily

plot_measure(dataset = merged, 
             measures ="TGFBR1", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="TGFBR2", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="TGFBR3",
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="TGFB1", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="BAMBI", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()


#NODAL/ActivinA

plot_measure(dataset = merged, 
             measures ="ACVR1B", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="ACVR2A", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="ACVR2B", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="FOXH1", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="SMAD2", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="SMAD3", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="SMAD4", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="NODAL", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="FURIN", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="TDGF1", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="LEFTY1", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="LEFTY2", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="CER1", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()


#Notch Genes


plot_measure(dataset = merged, 
             measures ="NOTCH1", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="NOTCH2", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="NOTCH3", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="DLL3", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged, 
             measures ="HES1", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()


plot_measure(dataset = merged, 
             measures ="RBPJ", 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
