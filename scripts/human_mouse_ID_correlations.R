library(Seurat)
library(tidyverse)

human <- readRDS("D:/realigning_for_int/signalling_objs/human/merged_integrated.RDS")
mouse <- readRDS("D:/realigning_for_int/signalling_objs/mouse/merged_mouse_hgncsymbols.rds")

human<-subset(human, idents='EPI')
mouse<-subset(mouse, idents='EPI')
DefaultAssay(human)<-'SCT'
DefaultAssay(mouse)<-'SCT'
human_matrix<-GetAssayData(human)
mouse_matrix<-GetAssayData(mouse)
human_meta<-as.data.frame(human@meta.data)
mouse_meta<-as.data.frame(mouse@meta.data)
rm(human)
rm(mouse)
gc()

stage='blastocyst'
human_blastocyst<-as.matrix(human_matrix[,human_meta$stage==stage])
human_gene<-as.numeric(human_blastocyst["ID2",])
human_correlations_blast<-apply(human_blastocyst,1,function(x){cor(human_gene,x)})
human_correlations_pval_blast<-apply(human_blastocyst,1,function(x){cor.test(human_gene,x)$p.value})
human_correlations_qval_blast=p.adjust(human_correlations_pval_blast,method = "BH")
human_significant_cor<-as.data.frame(human_correlations_qval_blast<0.05)
human_significant_cor<-subset(human_significant_cor, human_significant_cor[,1]==TRUE)
human_pearson_cors<-as.data.frame(human_correlations_blast)
human_pearson_cors<-as.data.frame(human_pearson_cors[rownames(human_significant_cor),])
rownames(human_pearson_cors)<-rownames(human_significant_cor)
colnames(human_pearson_cors)<-'blastocyst'


human_final<-data.frame(row.names=(unique(rownames(human_pearson_cors))))
human_final<-merge(human_final, human_pearson_cors, by='row.names', all.x=T)
rownames(human_final)<-human_final$Row.names
rm(human_pearson_cors)
gc()


stage='blastocyst'
mouse_blastocyst<-as.matrix(mouse_matrix[,mouse_meta$Stage==stage])
mouse_gene<-as.numeric(mouse_blastocyst["ID2",])
mouse_correlations_blast<-apply(mouse_blastocyst,1,function(x){cor(mouse_gene,x)})
mouse_correlations_pval_blast<-apply(mouse_blastocyst,1,function(x){cor.test(mouse_gene,x)$p.value})
mouse_correlations_qval_blast=p.adjust(mouse_correlations_pval_blast,method = "BH")
mouse_significant_cor<-as.data.frame(mouse_correlations_qval_blast<0.05)
mouse_significant_cor<-subset(mouse_significant_cor, mouse_significant_cor[,1]==TRUE)
mouse_pearson_cors<-as.data.frame(mouse_correlations_blast)
mouse_pearson_cors<-as.data.frame(mouse_pearson_cors[rownames(mouse_significant_cor),])
rownames(mouse_pearson_cors)<-rownames(mouse_significant_cor)
colnames(mouse_pearson_cors)<-'blastocyst'


mouse_final<-data.frame(row.names=(unique(rownames(mouse_pearson_cors))))
mouse_final<-merge(mouse_final, mouse_pearson_cors, by='row.names', all.x=T)
rownames(mouse_final)<-mouse_final$Row.names
rm(mouse_pearson_cors)

write.csv(human_final, 'human_ID2_sigcorrelations.csv')
write.csv(mouse_final, 'mouse_ID2_sigcorrelations.csv')

human_naive=c('SOX15', 'DNMT3L', 'KLF17', 'TFCP2L1', 'KLF4', 'ZFP42', 'PRDM14', 'TBX3','TFAP2C', 'DPPA3')
mouse_naive=c('NANOG', 'KLF2', 'KLF4', 'KLF17', 'PRDM14', 'DNMT3L', 'TBX3', 'DPPA3', 'SOX15', 'ZFP42')


human_naive_blast_dat=round(human_correlations_blast[human_naive],digits=2)
human_naive_blast_dat=human_naive_blast_dat[order(human_naive_blast_dat,decreasing=TRUE)]

mouse_naive_blast_dat=round(mouse_correlations_blast[mouse_naive],digits=2)
mouse_naive_blast_dat=mouse_naive_blast_dat[order(mouse_naive_blast_dat,decreasing=TRUE)]


pvals_vec_human_blast=human_correlations_qval_blast[c(names(human_naive_blast_dat))]
pvals_vec_mouse_blast=mouse_correlations_qval_blast[c(names(mouse_naive_blast_dat))]

human_pval=data.frame(human_blast=pvals_vec_human_blast)
write.table(human_pval,"pval_ID2_correlation_human.txt",quote=F,sep="\t")

mouse_pval=data.frame(mouse_blast=pvals_vec_mouse_blast)
write.table(mouse_pval,"pval_ID2_correlation_mouse.txt",quote=F,sep="\t")

#### Human Plots ###
#naive-blast
human_blast_dat=human_naive_blast_dat
cols=colorRampPalette(c("Blue","White","Red"))(201)
names(cols)=round(seq(-1,1,0.01),digits=2)
col_dat=cols[as.character(human_blast_dat)]
col_dat[is.na(col_dat)]="grey"
  
val_labels=human_blast_dat
val_labels[is.na(val_labels)]="N/A"
pdf(paste0("human_blast_ID2_correlations_naive.pdf"),useDingbats = F,width=14,height=3.5)
plot(0,type='n',axes=FALSE,ann=FALSE,ylim=c(0,2.1),xlim=c(-0.25,length(human_blast_dat)+0.2))
rect(xleft=seq(0,length(human_blast_dat)-1,1),xright=seq(1,length(human_blast_dat),1),ybottom=0,ytop=1,col=col_dat)
text(x=seq(0.35,length(human_blast_dat)-0.65,1),y=1.05,labels=names(human_blast_dat),srt=45,pos=4)
text(x=seq(0.5,length(human_blast_dat)-0.5,1),y=0.5,labels=val_labels)
text(x=-0.5,y=0.5,labels="ID2")
dev.off()

### Mouse Plots ###
#naive-blast
mouse_blast_dat=mouse_naive_blast_dat
cols=colorRampPalette(c("Blue","White","Red"))(201)
names(cols)=round(seq(-1,1,0.01),digits=2)
col_dat=cols[as.character(mouse_blast_dat)]
col_dat[is.na(col_dat)]="grey"
  
val_labels=mouse_blast_dat
val_labels[is.na(val_labels)]="N/A"
pdf(paste0("mouse_blast_ID2_correlations_naive.pdf"),useDingbats = F,width=14,height=3.5)
plot(0,type='n',axes=FALSE,ann=FALSE,ylim=c(0,2.1),xlim=c(-0.25,length(mouse_blast_dat)+0.2))
rect(xleft=seq(0,length(mouse_blast_dat)-1,1),xright=seq(1,length(mouse_blast_dat),1),ybottom=0,ytop=1,col=col_dat)
text(x=seq(0.35,length(mouse_blast_dat)-0.65,1),y=1.05,labels=names(mouse_blast_dat),srt=45,pos=4)
text(x=seq(0.5,length(mouse_blast_dat)-0.5,1),y=0.5,labels=val_labels)
text(x=-0.5,y=0.5,labels="ID2")
dev.off()



human_gene<-as.numeric(human_blastocyst["ID1",])
human_correlations_blast<-apply(human_blastocyst,1,function(x){cor(human_gene,x)})
human_correlations_pval_blast<-apply(human_blastocyst,1,function(x){cor.test(human_gene,x)$p.value})
human_correlations_qval_blast=p.adjust(human_correlations_pval_blast,method = "BH")
human_significant_cor<-as.data.frame(human_correlations_qval_blast<0.05)
human_significant_cor<-subset(human_significant_cor, human_significant_cor[,1]==TRUE)
human_pearson_cors<-as.data.frame(human_correlations_blast)
human_pearson_cors<-as.data.frame(human_pearson_cors[rownames(human_significant_cor),])
rownames(human_pearson_cors)<-rownames(human_significant_cor)
colnames(human_pearson_cors)<-'blastocyst'


human_final<-data.frame(row.names=(unique(rownames(human_pearson_cors))))
human_final<-merge(human_final, human_pearson_cors, by='row.names', all.x=T)
rownames(human_final)<-human_final$Row.names
rm(human_pearson_cors)
gc()


mouse_gene<-as.numeric(mouse_blastocyst["ID1",])
mouse_correlations_blast<-apply(mouse_blastocyst,1,function(x){cor(mouse_gene,x)})
mouse_correlations_pval_blast<-apply(mouse_blastocyst,1,function(x){cor.test(mouse_gene,x)$p.value})
mouse_correlations_qval_blast=p.adjust(mouse_correlations_pval_blast,method = "BH")
mouse_significant_cor<-as.data.frame(mouse_correlations_qval_blast<0.05)
mouse_significant_cor<-subset(mouse_significant_cor, mouse_significant_cor[,1]==TRUE)
mouse_pearson_cors<-as.data.frame(mouse_correlations_blast)
mouse_pearson_cors<-as.data.frame(mouse_pearson_cors[rownames(mouse_significant_cor),])
rownames(mouse_pearson_cors)<-rownames(mouse_significant_cor)
colnames(mouse_pearson_cors)<-'blastocyst'


mouse_final<-data.frame(row.names=(unique(rownames(mouse_pearson_cors))))
mouse_final<-merge(mouse_final, mouse_pearson_cors, by='row.names', all.x=T)
rownames(mouse_final)<-mouse_final$Row.names
rm(mouse_pearson_cors)

write.csv(human_final, 'human_ID1_sigcorrelations.csv')
write.csv(mouse_final, 'mouse_ID1_sigcorrelations.csv')

human_naive=c('SOX15', 'DNMT3L', 'KLF17', 'TFCP2L1', 'KLF4', 'ZFP42', 'PRDM14', 'TBX3', 'TFAP2C', 'DPPA3')
mouse_naive=c('NANOG', 'KLF2', 'KLF4', 'KLF17', 'PRDM14', 'DNMT3L', 'TBX3', 'DPPA3', 'SOX15', 'ZFP42')


human_naive_blast_dat=round(human_correlations_blast[human_naive],digits=2)
human_naive_blast_dat=human_naive_blast_dat[order(human_naive_blast_dat,decreasing=TRUE)]

mouse_naive_blast_dat=round(mouse_correlations_blast[mouse_naive],digits=2)
mouse_naive_blast_dat=mouse_naive_blast_dat[order(mouse_naive_blast_dat,decreasing=TRUE)]


pvals_vec_human_blast=human_correlations_qval_blast[c(names(human_naive_blast_dat))]
pvals_vec_mouse_blast=mouse_correlations_qval_blast[c(names(mouse_naive_blast_dat))]

human_pval=data.frame(human_blast=pvals_vec_human_blast)
write.table(human_pval,"pval_ID1_correlation_human.txt",quote=F,sep="\t")

mouse_pval=data.frame(mouse_blast=pvals_vec_mouse_blast)
write.table(mouse_pval,"pval_ID1_correlation_mouse.txt",quote=F,sep="\t")

#### Human Plots ###
#naive-blast
human_blast_dat=human_naive_blast_dat
cols=colorRampPalette(c("Blue","White","Red"))(201)
names(cols)=round(seq(-1,1,0.01),digits=2)
col_dat=cols[as.character(human_blast_dat)]
col_dat[is.na(col_dat)]="grey"
  
val_labels=human_blast_dat
val_labels[is.na(val_labels)]="N/A"
pdf(paste0("human_blast_ID1_correlations_naive.pdf"),useDingbats = F,width=14,height=3.5)
plot(0,type='n',axes=FALSE,ann=FALSE,ylim=c(0,2.1),xlim=c(-0.25,length(human_blast_dat)+0.2))
rect(xleft=seq(0,length(human_blast_dat)-1,1),xright=seq(1,length(human_blast_dat),1),ybottom=0,ytop=1,col=col_dat)
text(x=seq(0.35,length(human_blast_dat)-0.65,1),y=1.05,labels=names(human_blast_dat),srt=45,pos=4)
text(x=seq(0.5,length(human_blast_dat)-0.5,1),y=0.5,labels=val_labels)
text(x=-0.5,y=0.5,labels="ID1")
dev.off()

### Mouse Plots ###
#naive-blast
mouse_blast_dat=mouse_naive_blast_dat
cols=colorRampPalette(c("Blue","White","Red"))(201)
names(cols)=round(seq(-1,1,0.01),digits=2)
col_dat=cols[as.character(mouse_blast_dat)]
col_dat[is.na(col_dat)]="grey"
  
val_labels=mouse_blast_dat
val_labels[is.na(val_labels)]="N/A"
pdf(paste0("mouse_blast_ID1_correlations_naive.pdf"),useDingbats = F,width=14,height=3.5)
plot(0,type='n',axes=FALSE,ann=FALSE,ylim=c(0,2.1),xlim=c(-0.25,length(mouse_blast_dat)+0.2))
rect(xleft=seq(0,length(mouse_blast_dat)-1,1),xright=seq(1,length(mouse_blast_dat),1),ybottom=0,ytop=1,col=col_dat)
text(x=seq(0.35,length(mouse_blast_dat)-0.65,1),y=1.05,labels=names(mouse_blast_dat),srt=45,pos=4)
text(x=seq(0.5,length(mouse_blast_dat)-0.5,1),y=0.5,labels=val_labels)
text(x=-0.5,y=0.5,labels="ID1")
dev.off()




human_gene<-as.numeric(human_blastocyst["ID3",])
human_correlations_blast<-apply(human_blastocyst,1,function(x){cor(human_gene,x)})
human_correlations_pval_blast<-apply(human_blastocyst,1,function(x){cor.test(human_gene,x)$p.value})
human_correlations_qval_blast=p.adjust(human_correlations_pval_blast,method = "BH")
human_significant_cor<-as.data.frame(human_correlations_qval_blast<0.05)
human_significant_cor<-subset(human_significant_cor, human_significant_cor[,1]==TRUE)
human_pearson_cors<-as.data.frame(human_correlations_blast)
human_pearson_cors<-as.data.frame(human_pearson_cors[rownames(human_significant_cor),])
rownames(human_pearson_cors)<-rownames(human_significant_cor)
colnames(human_pearson_cors)<-'blastocyst'


human_final<-data.frame(row.names=(unique(rownames(human_pearson_cors))))
human_final<-merge(human_final, human_pearson_cors, by='row.names', all.x=T)
rownames(human_final)<-human_final$Row.names
rm(human_pearson_cors)
gc()


mouse_gene<-as.numeric(mouse_blastocyst["ID3",])
mouse_correlations_blast<-apply(mouse_blastocyst,1,function(x){cor(mouse_gene,x)})
mouse_correlations_pval_blast<-apply(mouse_blastocyst,1,function(x){cor.test(mouse_gene,x)$p.value})
mouse_correlations_qval_blast=p.adjust(mouse_correlations_pval_blast,method = "BH")
mouse_significant_cor<-as.data.frame(mouse_correlations_qval_blast<0.05)
mouse_significant_cor<-subset(mouse_significant_cor, mouse_significant_cor[,1]==TRUE)
mouse_pearson_cors<-as.data.frame(mouse_correlations_blast)
mouse_pearson_cors<-as.data.frame(mouse_pearson_cors[rownames(mouse_significant_cor),])
rownames(mouse_pearson_cors)<-rownames(mouse_significant_cor)
colnames(mouse_pearson_cors)<-'blastocyst'


mouse_final<-data.frame(row.names=(unique(rownames(mouse_pearson_cors))))
mouse_final<-merge(mouse_final, mouse_pearson_cors, by='row.names', all.x=T)
rownames(mouse_final)<-mouse_final$Row.names
rm(mouse_pearson_cors)

write.csv(human_final, 'human_ID3_sigcorrelations.csv')
write.csv(mouse_final, 'mouse_ID3_sigcorrelations.csv')

human_naive=c('SOX15', 'DNMT3L', 'KLF17', 'TFCP2L1', 'KLF4', 'ZFP42', 'PRDM14', 'TBX3', 'TFAP2C', 'DPPA3')
mouse_naive=c('NANOG', 'KLF2', 'KLF4', 'KLF17', 'PRDM14', 'DNMT3L', 'TBX3', 'DPPA3', 'SOX15', 'ZFP42')


human_naive_blast_dat=round(human_correlations_blast[human_naive],digits=2)
human_naive_blast_dat=human_naive_blast_dat[order(human_naive_blast_dat,decreasing=TRUE)]

mouse_naive_blast_dat=round(mouse_correlations_blast[mouse_naive],digits=2)
mouse_naive_blast_dat=mouse_naive_blast_dat[order(mouse_naive_blast_dat,decreasing=TRUE)]


pvals_vec_human_blast=human_correlations_qval_blast[c(names(human_naive_blast_dat))]
pvals_vec_mouse_blast=mouse_correlations_qval_blast[c(names(mouse_naive_blast_dat))]

human_pval=data.frame(human_blast=pvals_vec_human_blast)
write.table(human_pval,"pval_ID3_correlation_human.txt",quote=F,sep="\t")

mouse_pval=data.frame(mouse_blast=pvals_vec_mouse_blast)
write.table(mouse_pval,"pval_ID3_correlation_mouse.txt",quote=F,sep="\t")

#### Human Plots ###
#naive-blast
human_blast_dat=human_naive_blast_dat
cols=colorRampPalette(c("Blue","White","Red"))(201)
names(cols)=round(seq(-1,1,0.01),digits=2)
col_dat=cols[as.character(human_blast_dat)]
col_dat[is.na(col_dat)]="grey"
  
val_labels=human_blast_dat
val_labels[is.na(val_labels)]="N/A"
pdf(paste0("human_blast_ID3_correlations_naive.pdf"),useDingbats = F,width=14,height=3.5)
plot(0,type='n',axes=FALSE,ann=FALSE,ylim=c(0,2.1),xlim=c(-0.25,length(human_blast_dat)+0.2))
rect(xleft=seq(0,length(human_blast_dat)-1,1),xright=seq(1,length(human_blast_dat),1),ybottom=0,ytop=1,col=col_dat)
text(x=seq(0.35,length(human_blast_dat)-0.65,1),y=1.05,labels=names(human_blast_dat),srt=45,pos=4)
text(x=seq(0.5,length(human_blast_dat)-0.5,1),y=0.5,labels=val_labels)
text(x=-0.5,y=0.5,labels="ID3")
dev.off()

### Mouse Plots ###
#naive-blast
mouse_blast_dat=mouse_naive_blast_dat
cols=colorRampPalette(c("Blue","White","Red"))(201)
names(cols)=round(seq(-1,1,0.01),digits=2)
col_dat=cols[as.character(mouse_blast_dat)]
col_dat[is.na(col_dat)]="grey"
  
val_labels=mouse_blast_dat
val_labels[is.na(val_labels)]="N/A"
pdf(paste0("mouse_blast_ID3_correlations_naive.pdf"),useDingbats = F,width=14,height=3.5)
plot(0,type='n',axes=FALSE,ann=FALSE,ylim=c(0,2.1),xlim=c(-0.25,length(mouse_blast_dat)+0.2))
rect(xleft=seq(0,length(mouse_blast_dat)-1,1),xright=seq(1,length(mouse_blast_dat),1),ybottom=0,ytop=1,col=col_dat)
text(x=seq(0.35,length(mouse_blast_dat)-0.65,1),y=1.05,labels=names(mouse_blast_dat),srt=45,pos=4)
text(x=seq(0.5,length(mouse_blast_dat)-0.5,1),y=0.5,labels=val_labels)
text(x=-0.5,y=0.5,labels="ID3")
dev.off()

