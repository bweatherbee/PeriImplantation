library(Seurat)


merged <- readRDS("D:/realigning_for_int/FOR 2022 SIGNALING PAPER/human/merged_integrated.RDS")
DefaultAssay(merged)<-"SCT"

blastocyst<-subset(merged, subset=stage=="blastocyst")
Idents(blastocyst)<-blastocyst$new_cell_type
blastocyst<-subset(blastocyst, idents = c("EPI", "HYPO", "TE"))

blast_matrix<-blastocyst@assays[["SCT"]]@data
blast_matrix<-as.matrix(blast_matrix)
blast_meta<-cbind(blastocyst@meta.data[,'new_cell_type', drop=F])
write.table(blast_matrix, './blastocyst/blast_matrix.txt', sep='\t', quote=F)
write.table(blast_meta, './blastocyst/blast_meta.txt', sep='\t', quote=F)

peri<-subset(merged, subset=stage=="peri")
peri_matrix<-as.matrix(peri@assays[["SCT"]]@data)
peri_meta<-cbind(peri@meta.data[,'new_cell_type', drop=F])
write.table(peri_matrix, './peri/peri_matrix.txt', sep='\t', quote=F)
write.table(peri_meta, './peri/peri_meta.txt', sep='\t', quote=F)

e_post<-subset(merged, subset=stage=="early post")
epost_matrix<-as.matrix(e_post@assays[["SCT"]]@data)
epost_meta<-cbind(e_post@meta.data[,'new_cell_type', drop=F])
write.table(epost_matrix, './e_post/epost_matrix.txt', sep='\t', quote=F)
write.table(epost_meta, './e_post/epost_meta.txt', sep='\t', quote=F)

l_post<-subset(merged, subset=stage=="late post")
lpost_matrix<-as.matrix(l_post@assays[["SCT"]]@data)
lpost_meta<-cbind(l_post@meta.data[,'new_cell_type', drop=F])
write.table(lpost_matrix, './l_post/lpost_matrix.txt', sep='\t', quote=F)
write.table(lpost_meta, './l_post/lpost_meta.txt', sep='\t', quote=F)



#In text editor, add column name "cell" to meta file
#in ubuntu window, install cellphonedb and R4.0
#in ubuntu window, change to directory with counts & run
  #cellphonedb method statistical_analysis meta_file.txt matrix_file.txt --iterations=10 --threads=2 --counts-data=gene_name
  #cellphonedb plot dot_plot --means-path=./out/means.txt --pvalues-path=./out/pvalues.txt --rows rows_to_plot.txt
