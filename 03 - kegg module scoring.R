library(Seurat)
library(biomaRt)
library(org.Hs.eg.db)
#BiocManager::install("EnrichmentBrowser")
library(EnrichmentBrowser)

merged <- readRDS("D:/realigning_for_int/FOR 2022 SIGNALING PAPER/human/merged_integrated.RDS")

kegg<-getGenesets(org='hsa', db="kegg", gene.id.type="ENTREZID")


ensembl=useMart("ensembl")
datasets<-listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl", mart=ensembl)
MAPK<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
             filters="entrezgene_id",
             values=kegg$hsa04010_MAPK_signaling_pathway, mart=ensembl)

mTOR<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="entrezgene_id",
            values=kegg$hsa04150_mTOR_signaling_pathway, mart=ensembl)

PI3K_AKT<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="entrezgene_id",
            values=kegg$`hsa04151_PI3K-Akt_signaling_pathway`, mart=ensembl)
WNT<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="entrezgene_id",
            values=kegg$hsa04310_Wnt_signaling_pathway, mart=ensembl)
Notch<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="entrezgene_id",
            values=kegg$hsa04330_Notch_signaling_pathway, mart=ensembl)
Hedgehog<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="entrezgene_id",
            values=kegg$hsa04340_Hedgehog_signaling_pathway, mart=ensembl)
TGFB<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="entrezgene_id",
            values=kegg$`hsa04350_TGF-beta_signaling_pathway`, mart=ensembl)
Hippo1<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="entrezgene_id",
            values=kegg$hsa04390_Hippo_signaling_pathway, mart=ensembl)
Hippo2<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="entrezgene_id",
            values=kegg$hsa04392_Hippo_signaling_pathway, mart=ensembl)
Jak_Stat<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters="entrezgene_id",
            values=kegg$`hsa04630_JAK-STAT_signaling_pathway`, mart=ensembl)

MAPK<-list(unique(MAPK$external_gene_name))
mTOR<-list(unique(mTOR$external_gene_name))
PI3K_AKT<-list(unique(PI3K_AKT$external_gene_name))
WNT<-list(unique(WNT$external_gene_name))
Notch<-list(unique(Notch$external_gene_name))
Hedgehog<-list(unique(Hedgehog$external_gene_name))
TGFB<-list(unique(TGFB$external_gene_name))
Hippo1<-list(unique(Hippo1$external_gene_name))
Hippo2<-list(unique(Hippo2$external_gene_name))
Jak_Stat<-list(unique(Jak_Stat$external_gene_name))


DefaultAssay(merged)<-"SCT"
merged_scored<-AddModuleScore(merged, features=MAPK, assay="UMI_qumi", name="MAPK")
merged_scored$MAPK<-merged_scored$MAPK1
merged_scored$MAPK1<-NULL

merged_scored<-AddModuleScore(merged_scored, features=mTOR, assay="UMI_qumi", name="mTOR")
merged_scored$mTOR<-merged_scored$mTOR1
merged_scored$mTOR1<-NULL

merged_scored<-AddModuleScore(merged_scored, features=PI3K_AKT, assay="UMI_qumi", name="PI3K_AKT")
merged_scored$PI3K_AKT<-merged_scored$PI3K_AKT1
merged_scored$PI3K_AKT1<-NULL

merged_scored<-AddModuleScore(merged_scored, features=WNT, assay="UMI_qumi", name="WNT")
merged_scored$WNT<-merged_scored$WNT1
merged_scored$WNT1<-NULL

merged_scored<-AddModuleScore(merged_scored, features=Notch, assay="UMI_qumi", name="Notch")
merged_scored$Notch<-merged_scored$Notch1
merged_scored$Notch1<-NULL

merged_scored<-AddModuleScore(merged_scored, features=Hedgehog, assay="UMI_qumi", name="Hedgehog")
merged_scored$Hedgehog<-merged_scored$Hedgehog1
merged_scored$Hedgehog1<-NULL

merged_scored<-AddModuleScore(merged_scored, features=TGFB, assay="UMI_qumi", name="TGFB")
merged_scored$TGFB<-merged_scored$TGFB1
merged_scored$TGFB1<-NULL

merged_scored<-AddModuleScore(merged_scored, features=Hippo1, assay="UMI_qumi", name="Hippo1")
merged_scored$Hippo1<-merged_scored$Hippo11
merged_scored$Hippo11<-NULL

merged_scored<-AddModuleScore(merged_scored, features=Hippo2, assay="UMI_qumi", name="Hippo2")
merged_scored$Hippo2<-merged_scored$Hippo21
merged_scored$Hippo21<-NULL

merged_scored<-AddModuleScore(merged_scored, features=Jak_Stat, assay="UMI_qumi", name="Jak_Stat")
merged_scored$Jak_Stat<-merged_scored$Jak_Stat1
merged_scored$Jak_Stat1<-NULL



library(Scillus)
library(RColorBrewer)

merged_scored@assays$RNA <- merged_scored@assays$UMI_qumi

plot_measure(dataset = merged_scored, 
             measures = c("MAPK", "mTOR", "PI3K_AKT", "WNT", "Notch", "Hedgehog", "TGFB", "Hippo1", "Hippo2", "Jak_Stat"), 
             group_by = "new_cell_type",
             pal_setup = "Dark2")

plot_measure(dataset = merged_scored, 
             measures = c("MAPK"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged_scored, 
             measures = c("mTOR"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged_scored, 
             measures = c("PI3K_AKT"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged_scored, 
             measures = c("WNT"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1",
             plot_type="combined")+theme_bw()

plot_measure(dataset = merged_scored, 
             measures = c("Notch"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged_scored, 
             measures = c("Hedgehog"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged_scored, 
             measures = c("TGFB"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged_scored, 
             measures = c("Hippo1"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged_scored, 
             measures = c("Hippo2"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure(dataset = merged_scored, 
             measures = c("Jak_Stat"), 
             group_by = "stage",
             split_by = "new_cell_type",
             pal_setup = "Pastel1")+theme_bw()

plot_measure_dim(dataset = merged_scored, 
                 measures = c("MAPK", "mTOR", "PI3K_AKT", "WNT", "Notch", "Hedgehog", "TGFB", "Hippo1", "Hippo2", "Jak_Stat"))

FeaturePlot(merged_scored, features=c("MAPK", "mTOR", "PI3K_AKT", "WNT", "Notch", "Hedgehog", "TGFB", "Hippo1", "Hippo2", "Jak_Stat"),
            min.cutoff = 'q5', max.cutoff = 'q98', order=TRUE, cols=c("grey88", "green3"), pt.size=2)


metadata<-merged_scored@meta.data
signaling_matrix<-metadata[,24:33]
Scored<-CreateSeuratObject(counts=t(signaling_matrix),assay='signaling', meta.data = merged@meta.data)

Scored<-ScaleData(Scored)
plot_heatmap(dataset=Scored, markers=c("MAPK", "mTOR", "PI3K-AKT", "WNT", "Notch", "Hedgehog", "TGFB", "Hippo1", "Hippo2", "Jak-Stat"),
             sort_var=c("new_cell_type", "Age"),
             anno_var=c("new_cell_type", "Age"),
             anno_colors = c("Set1", "Set2", "Set3"),
             hm_limit = c(-1,0,1))
Idents(Scored)<-Scored$new_cell_type
avg_activity<-as.data.frame(AverageExpression(Scored, assays="signaling", slot="counts", group.by = "new_cell_type", add.ident="Age"))
pheatmap::pheatmap(t(scale(t(avg_activity))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
avg_activity2<-as.data.frame(AverageExpression(Scored, assays="signaling", slot="counts", group.by="new_cell_type"))
pheatmap::pheatmap(t(scale(t(avg_activity2))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

avg_activity3<-as.data.frame(AverageExpression(Scored, assays="signaling", slot="counts", group.by="new_cell_type", add.ident="stage"))
pheatmap::pheatmap(t(scale(t(avg_activity3))), cluster_cols = FALSE,  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
