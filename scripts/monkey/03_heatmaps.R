library(Seurat)
library(Scillus)
merged <- readRDS("C:/Users/baile/Desktop/realigning_for_int/CynMonkey/merged_hgnc19_symbols.rds")
DefaultAssay(merged)<-"RNA"
merged<-ScaleData(merged, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features=rownames(merged))
saveRDS(merged, 'merged_hgnc19_symbols.rds')

BMP_components<-c("BMP1", "BMP2", "BMP3", "BMP4", "BMP5", "BMP6", "BMP7",
                  "BMP8A","BMP8B", "GDF2", "BMP10", "GDF11", "GDF7", "GDF6",
                  "SMAD1", "SMAD5", "SMAD9", "SMAD4", "SMAD6", "SMAD7",
                  "BMPR1A", "BMPR1B", "ACVR1", "BMPR2", "ACVR2A", "ACVR2B",
                  "ID1", "ID2", "ID3", "ID4")

NODAL_components<-c("NODAL", "GDF3",
            "CFC1B", "TDGF1", "ACVR1B", "ACVR2A", "ACVR2B", 
                    "LEFTY1", "LEFTY2", "CER1",
                    "SMAD2", "SMAD3", "SMAD4",
                    "FOXH1", "MIXL1")

WNT_components<-c("WNT1", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT9A", "WNT9B", "WNT10A", "WNT10B", "WNT11", "WNT16",
                  "FZD1", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "FZD10",
                  "LRP5", "LRP6", "LRP1", "KREMEN1", "KREMEN2", "ROR2", "CELSR1", "CELSR2", "CELSR3",
                  "DKK1", "DKK2", "DKK3", "DKK4", "SOSTDC1", "WIF1", "CER1", "SFRP1", "SFRP2", "FRZB", "SFRP4", "SFRP5",
                  "TCF7", "LEF1", "AXIN2")

Notch_components<-c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
                    "JAG1", "JAG2", "DLL1", "DLL3", "DLL4",
                    "DLK1", "DLK2")


levels(merged$new_cell_type)
my_levels<-c("EPI", "HYPO", "TE", "EXMC")
merged$new_cell_type <- factor(x = merged$new_cell_type, levels = my_levels)

plot_heatmap(dataset=merged, markers=BMP_components,
             sort_var=c("new_cell_type", "Age", "paper"),
             anno_var=c("new_cell_type", "Age", "paper"),
             anno_colors = c("Set1", "Set2", "Set3"),
             hm_limit = c(-1,0,1))

plot_heatmap(dataset=merged, markers=NODAL_components,
             sort_var=c("new_cell_type", "Age", "paper"),
             anno_var=c("new_cell_type", "Age", "paper"),
             anno_colors = c("Set1", "Set2", "Set3"),
             hm_limit = c(-1,0,1))

plot_heatmap(dataset=merged, markers=WNT_components,
             sort_var=c("new_cell_type", "Age", "paper"),
             anno_var=c("new_cell_type", "Age", "paper"),
             anno_colors = c("Set1", "Set2", "Set3"),
             hm_limit = c(-1,0,1))

plot_heatmap(dataset=merged, markers=Notch_components,
             sort_var=c("new_cell_type", "Age", "paper"),
             anno_var=c("new_cell_type", "Age", "paper"),
             anno_colors = c("Set1", "Set2", "Set3"),
             hm_limit = c(-1,0,1))

DimPlot(merged, group.by='new_cell_type')
DimPlot(merged, group.by='paper')
DimPlot(merged, group.by='Age')
DimPlot(merged, group.by='stage')


FeaturePlot(merged, features=c("SOX2", "NANOG", "POU5F1", "KLF17"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("SOX17", "GATA4", "GATA6", "PDGFRA"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("GATA3", "GATA2","TFAP2C", "CDX2"),order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("SDC1", "CGB1", "CGA", "ITGB3"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("CER1", "LEFTY1", "LEFTY2", "DKK1"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')


VlnPlot(merged, features=BMP_components, group.by = "new_cell_type")
VlnPlot(merged, features=NODAL_components, group.by = "new_cell_type")
VlnPlot(merged, features=WNT_components, group.by = "new_cell_type")
VlnPlot(merged, features=Notch_components, group.by = "new_cell_type")

### Vlns of selected Genes ###
DefaultAssay(merged)<-'SCT'
levels(merged$stage)
my_levels<-c("blastocyst", "peri", "early post", "late post", "gast")
merged$stage <- factor(x = merged$stage, levels = my_levels)


VlnPlot(merged, features=c('NODAL', 'FOXH1', 'TDGF1', 'PCSK6', 'LEFTY1', 'LEFTY2',
                           'BMPR1A', 'BMPR1B', 'BMPR2', 'ID1', 'ID2', 'ID3', 'ID4',
                           'NOTCH1', 'DLL1', 'DLL3', 'JAG1', 'HES1', 'RBPJ',
                           'FGF2', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3','FGFR4',
                           'WNT3', 'WNT5A', 'TCF7', 'TCF7L1', 'TCF7L2', 'LEF1'), group.by='new_cell_type',
        split.by = 'stage', 
        cols= RColorBrewer::brewer.pal(5, 'Pastel1'), pt.size=0, stack=T)

DefaultAssay(merged)<-'RNA'
merged<-NormalizeData(merged, normalization.method = 'LogNormalize')
VlnPlot(merged, features=c('NODAL', 'FOXH1', 'TDGF1', 'PCSK6', 'LEFTY1', 'LEFTY2',
                          'BMPR1A', 'BMPR1B', 'BMPR2', 'ID1', 'ID2', 'ID3', 'ID4',
                          'NOTCH1', 'DLL1', 'DLL3', 'JAG1', 'HES1', 'RBPJ',
                          'FGF2', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3','FGFR4',
                          'WNT3', 'WNT5A', 'TCF7', 'TCF7L1', 'TCF7L2', 'LEF1'), group.by='new_cell_type',
       split.by = 'stage', 
       cols= RColorBrewer::brewer.pal(5, 'Pastel1'), pt.size=0, stack=T)