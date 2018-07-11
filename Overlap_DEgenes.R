####Load necessities

library(robustbase)
library(dplyr)
require(ggrepel)
require(pheatmap)
require(VennDiagram)
library(base)
require(gplots)
require(VennDiagram)
require(xlsx)

source("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Required functions.R")

setwd("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D")

load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/pancreas.Robj")

load("significantgenes_islets.Robj")
islets <- genes
load("significantgenes_alphacells.Robj")
alphagenes <- genes
load("significantgenes_betacells.Robj")
betagenes <- genes

#### Overlap alpha and beta (sc)
title = "T2D vs Healthy: Alpha and beta signature genes"
venn <- venn.diagram(x=list(betagenes, alphagenes), "Overlap sc_ab.png", imagetype="png", category.names = c("Sc-Beta cells", "Sc-Alpha cells"), main=title, fill=c("#7570b3", "#e7298a"), cat.pos=c(-25,150))
overlap_ab <- calculate.overlap(list(alphagenes, betagenes))  
shared_sc <- overlap_ab$a3
shared_sc<-(unique(sub("(\\_|\\.).+","", shared_sc)))
unique_alpha <- setdiff(x=alphagenes, y=betagenes)
unique_alpha<-(unique(sub("(\\_|\\.).+","", unique_alpha)))
unique_beta <- setdiff(x=betagenes, y=alphagenes)
unique_beta<-(unique(sub("(\\_|\\.).+","", unique_beta)))


####Overlap bulk and sc
title = "T2D vs Healthy: Overlap bulk and sc-seq signature genes"
venn <- venn.diagram(x=list(betagenes, alphagenes, islets), "Overlap sc_islets.png", imagetype="png", category.names = c("Sc-Beta cells", "Sc-Alpha cells", "Islets"), main=title, fill=c("#7570b3", "#e7298a", "#e6ab02"), cat.pos=c(-40,150,0))
overlap_bulksc <- calculate.overlap(list(alphagenes,betagenes,islets))
shared_bulksc <- overlap_bulksc$a5
shared_sc_top10 <- Wx2[shared_bulksc,]
shared_sc_top10<-shared_bulksc_top10[order(abs(shared_bulksc_top10$lfc2),decreasing = T),]


shared_sc_top10 <- shared_bulksc_top10[1:10,]
sc <- c(alphagenes, betagenes)
unique_bulk <- setdiff(islets, sc)
unique_bulk<-(unique(sub("(\\_|\\.).+","", unique_bulk)))
missing_betab <- overlap_bulksc$a3
missing_alphab <- overlap_bulksc$a1

####Save 
write.csv(shared_sc,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/shared_sc.csv",quote=F,row.names=F)

write.csv(unique_alpha,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/unique_alphagenes.csv",quote=F,row.names=F)

write.csv(unique_beta,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/unique_betagenes.csv",quote=F,row.names=F)

write.csv(unique_bulk,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/unique_bulkgenes.csv",quote=F,row.names=F)

write.csv(missing_betab,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/betagenes_missedbulk.csv",quote=F,row.names=F)

write.csv(missing_alphab,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/alphagenes_missedbulk.csv",quote=F,row.names=F)

write.csv(shared_sc_top10, "Top 10 shared sc")