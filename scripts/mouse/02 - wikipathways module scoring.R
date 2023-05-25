library(Seurat)
library(biomaRt)
library(org.Hs.eg.db)
#BiocManager::install("rWikiPathways")
library(rWikiPathways)

merged <- readRDS("D:/realigning_for_int/FOR 2022 SIGNALING PAPER/mouse/merged_mouse_hgncsymbols.rds")
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
signaling_matrix<-merged_scored@meta.data[,29:55]
colnames(signaling_matrix)<-c("AMPK", "BMP1", "BMP2", "NODAL", "Activin", "TGFBR_C", "TGFB", "TGFBR",
                              "IGF1R", "Notch", "Hippo", "Hippo_Yap", "Hippo_reg", "MAPK", "MAPK4_6",
                              "FGFR1", "FGFR2", "FGFR3", "FGFR4", "PI3K_AKT", "mTOR", "Hh", "RA",
                              "WNT1", "WNT2", "WNT_bcat_ind", "WNT_TCF_dep")
merged_scored<-AddMetaData(merged, metadata=signaling_matrix)
saveRDS(signaling_matrix, 'wikipathways_modulescores.rds')


library(Scillus)
library(RColorBrewer)
library(ggplot2)
merged_scored@assays[['RNA']]<-merged_scored@assays$RNA_qumi
plot_measure(dataset = merged_scored, 
             measures = c("AMPK", "BMP1", "BMP2", "NODAL", "Activin", "TGFBR_C", "TGFB", "TGFBR",
                          "IGF1R", "Notch", "Hippo", "Hippo_Yap", "Hippo_reg", "MAPK", "MAPK4_6",
                          "FGFR1", "FGFR2", "FGFR3", "FGFR4", "PI3K_AKT", "mTOR", "Hh", "RA",
                          "WNT1", "WNT2", "WNT_bcat_ind", "WNT_TCF_dep"), 
             group_by = "new_cell_type",
             pal_setup = "Dark2")

plot_measure(dataset = merged_scored, 
             measures = c("AMPK"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("BMP1"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1") +theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("BMP2"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("NODAL"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("Activin"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("TGFBR_C"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("TGFB"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("TGFBR"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("IGF1R"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("Notch"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("Hippo"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("Hippo_Yap"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("Hippo_reg"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("MAPK"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("MAPK4_6"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("FGFR1"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()
plot_measure(dataset = merged_scored, 
             measures = c("FGFR2"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("FGFR3"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("FGFR4"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("PI3K_AKT"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("mTOR"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("Hh"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("RA"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("WNT1"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("WNT2"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("WNT_bcat_ind"), 
             group_by = "Stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")
plot_measure(dataset = merged_scored, 
             measures = c("WNT_TCF_dep"), 
             group_by = "Stage",
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
             sort_var=c("new_cell_type", "Stage", "Paper"),
             anno_var=c("new_cell_type", "Stage", "Paper"),
             anno_colors = c("Set1", "Set2", "Set3"),
             hm_limit = c(-1,0,1))

Idents(Scored)<-Scored$new_cell_type
avg_activity<-as.data.frame(AverageExpression(Scored, assays="signaling", slot="counts", group.by = "new_cell_type", add.ident="Age"))
pheatmap::pheatmap(t(scale(t(avg_activity))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
avg_activity2<-as.data.frame(AverageExpression(Scored, assays="signaling", slot="counts", group.by="new_cell_type"))
pheatmap::pheatmap(t(scale(t(avg_activity2))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))


metadata<-merged@meta.data

Idents(Scored)<-Scored$new_cell_type

avg_activity3<-as.data.frame(AverageExpression(Scored, assays="signaling", slot="counts", group.by="new_cell_type", add.ident="Stage"))

pheatmap::pheatmap(t(scale(t(avg_activity3))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

DefaultAssay(merged)<-'SCT'
levels(merged$Stage)

#Plots for comparative violins of selected genes
VlnPlot(merged, features=c('Nodal', 'Foxh1', 'Cfc1', 'Pcsk6', 'Lefty1', 'Lefty2',
                           'Bmpr1a', 'Bmpr1b', 'Bmpr2', 'Id1', 'Id2', 'Id3', 'Id4',
                           'Notch1', 'Dll1', 'Dll3', 'Jag1', 'Hes1', 'Rbpj',
                           'Fgf2', 'Fgf4', 'Fgfr1', 'Fgfr2', 'Fgfr3','Fgfr4',
                           'Wnt3', 'Wnt5a', 'Tcf7', 'Tcf7l1', 'Tcf7l2', 'Lef1'), group.by='new_cell_type',
        split.by = 'Stage', 
        cols= c("#FBB4AE", "#CCEBC5", "#DECBE4", "#FED9A6"), pt.size=0, stack=T)

DefaultAssay(merged)<-'RNA_qumi'
merged<-NormalizeData(merged, normalization.method = 'LogNormalize')
VlnPlot(merged, features=c('Nodal', 'Foxh1', 'Cfc1', 'Pcsk6', 'Lefty1', 'Lefty2',
                           'Bmpr1a', 'Bmpr1b', 'Bmpr2', 'Id1', 'Id2', 'Id3', 'Id4',
                           'Notch1', 'Dll1', 'Dll3', 'Jag1', 'Hes1', 'Rbpj',
                           'Fgf2', 'Fgf4', 'Fgfr1', 'Fgfr2', 'Fgfr3','Fgfr4',
                           'Wnt3', 'Wnt5a', 'Tcf7', 'Tcf7l1', 'Tcf7l2', 'Lef1'), group.by='new_cell_type',
        split.by = 'Stage', 
        cols= c("#FBB4AE", "#CCEBC5", "#DECBE4", "#FED9A6"), pt.size=0, stack=T)

