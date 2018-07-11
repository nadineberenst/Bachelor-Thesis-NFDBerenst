####Load all necessities

setwd("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/ValidationID")

load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/pancreas.Robj")


source("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Required functions.R")

library(dplyr)
library(ggrepel)
library(robustbase)
library(Seurat)
library(base)

####Find all significant genes 
healthy <- pancreas@data.info[(pancreas@data.info$Condition == "Healthy"),]
healthy <- grep("P546", rownames(healthy), invert=T, value=T) 
healthy <- pancreas@data.info[healthy,]  

donors <-  c("D1", "D2", "D3", "D4", "D6", "D7", "D8")

for (i in 1:length(donors)){
  
  fname <- paste(donors[i], "test", sep="-")
  dir.create(fname, showWarnings = FALSE)
    train_info <- healthy[grep(donors[i], healthy$Donor, invert=T),]
    train_info <- na.omit(train_info)
    
a <- rownames(train_info[(train_info$CT_TSNE == "Beta"),]); a2= "Beta"
b <- rownames(train_info[(train_info$CT_TSNE == "Alpha"),]); b2 = "Alpha"
  
A<-data.matrix(pancreas@data[,a]); dim(A)
B<-data.matrix(pancreas@data[,b]); dim(B)

Wx<- sapply(1:nrow(A), function(i){
  wilcox.test(as.numeric(A[i,]), as.numeric(B[i,]), exact=FALSE)$p.value
}) 

Wx2<-data.frame(row.names = rownames(A),p_val=Wx)

Wx2<- cbind(Wx2,Padj=p.adjust(Wx2$p_val, "BH"),
            lfc1=log2(rowMedians(as.matrix(pancreas@data[rownames(Wx2),a])))-
              log2(rowMedians(as.matrix(pancreas@data[rownames(Wx2),b]))),
            lfc2=log2(rowMeans(as.matrix(pancreas@data[rownames(Wx2),a])))-
              log2(rowMeans(as.matrix(pancreas@data[rownames(Wx2),b]))),
            rowMedians(as.matrix(pancreas@data[rownames(Wx2),a])),
            rowMedians(as.matrix(pancreas@data[rownames(Wx2),b])),
            rowMeans(as.matrix(pancreas@data[rownames(Wx2),a])),
            rowMeans(as.matrix(pancreas@data[rownames(Wx2),b]))
)
colnames(Wx2)[c(5,6,7,8)]<-c(paste("Med",a2,sep="_"),paste("Med",b2,sep="_"),paste("Av",a2,sep="_"),paste("Av",b2,sep="_"))

save(Wx2, file= paste(fname, "Wx2", sep="/"))

genes <- rownames(Wx2[Wx2$Padj<0.01 & abs(Wx2$lfc2)>1,])
 
save(genes, file = paste(fname, "DG", sep="/"))   

train_data <- pancreas@data[,rownames(train_info)]
pancdata<- train_data / rowMeans(train_data)

A<-rowMeans(pancdata[genes,b])
B<-rowMeans(pancdata[genes,a])
C<-pancdata[genes,c(b,a)]
AB<-A-B
AB<-na.omit(AB)
CB<-sweep(C,1,B,"-");
CB<-na.omit(CB)
t1<-sweep(CB,1,AB,"*")
t1<-colSums(t1)
t2<-sqrt(sum(AB^2))
score<-t1/(t2^2)
score<-score[order(score, decreasing=F)]
score_all <- score


save(score_all, file=paste(fname, "score_all_train.Robj", sep="/"))
save(pancdata, file=paste(fname, "pancdata_train.Robj", sep="/"))
save(A, file=paste(fname, "score DG Alpha_train.Robj", sep="/"))
save(B, file=paste(fname, "score DG Beta_train.Robj", sep="/"))

Amean<-mean(score_all[b])
Asd<-sd(score_all[b])
Aco<-Amean-2*Asd
Aco_plus <- Amean+2*Asd

Bmean<-mean(score_all[a])
Bsd<-sd(score_all[a])
Bco<-Bmean+2*Bsd
Bco_min <- Bmean-2*Bsd

save(Amean, file=paste(fname, "score Alpha mean_all.Robj", sep="/"))
save(Bmean, file=paste(fname, "score Beta mean_all.Robj", sep="/"))
save(Aco, file=paste(fname, "score Aco_all.Robj", sep="/"))
save(Aco_plus, file=paste(fname, "score Aco+_all.Robj", sep="/"))
save(Bco, file=paste(fname, "score Bco_all.Robj", sep="/"))
save(Bco_min, file=paste(fname, "score Bco-_all.Robj", sep="/"))


}

####Calculate ID-scores for overlapping DGs
healthy <- pancreas@data.info[(pancreas@data.info$Condition == "Healthy"),]
healthy <- grep("P546", rownames(healthy), invert=T, value=T) 
healthy <- pancreas@data.info[healthy,]  
donors <-  c("D1", "D2", "D3", "D4", "D6", "D7", "D8")

for (i in 1:length(donors)){
  
  fname <- paste(donors[i], "test", sep="-")
  dir.create(fname, showWarnings = FALSE)
  train_info <- healthy[grep(donors[i], healthy$Donor, invert=T),]
  train_info <- na.omit(train_info)
  
  a <- rownames(train_info[(train_info$CT_TSNE == "Beta"),]); a2= "Beta"
  b <- rownames(train_info[(train_info$CT_TSNE == "Alpha"),]); b2 = "Alpha"
  
  A<-data.matrix(pancreas@data[,a]); dim(A)
  B<-data.matrix(pancreas@data[,b]); dim(B)
  
  load(paste(fname,"DG_overlap", sep="/"))
  
  genes <- overlap_genes
 
  train_data <- pancreas@data[,rownames(train_info)]
  pancdata<- train_data / rowMeans(train_data)
  
  A<-rowMeans(pancdata[genes,b])
  B<-rowMeans(pancdata[genes,a])
  C<-pancdata[genes,c(b,a)]
  AB<-A-B
  AB<-na.omit(AB)
  CB<-sweep(C,1,B,"-");
  CB<-na.omit(CB)
  t1<-sweep(CB,1,AB,"*")
  t1<-colSums(t1)
  t2<-sqrt(sum(AB^2))
  score<-t1/(t2^2)
  score<-score[order(score, decreasing=F)]
  score_all <- score
  
  
  save(score_all, file=paste(fname, "score_all_train_overlap.Robj", sep="/"))
  save(pancdata, file=paste(fname, "pancdata_train_overlap.Robj", sep="/"))
  save(A, file=paste(fname, "score DG Alpha_train_overlap.Robj", sep="/"))
  save(B, file=paste(fname, "score DG Beta_train_overlap.Robj", sep="/"))
  
  Amean<-mean(score_all[b])
  Asd<-sd(score_all[b])
  Aco<-Amean-2*Asd
  Aco_plus <- Amean+2*Asd
  
  Bmean<-mean(score_all[a])
  Bsd<-sd(score_all[a])
  Bco<-Bmean+2*Bsd
  Bco_min <- Bmean-2*Bsd
  
  save(Amean, file=paste(fname, "score Alpha mean_train_overlap.Robj", sep="/"))
  save(Bmean, file=paste(fname, "score Beta mean_train_overlap.Robj", sep="/"))
  save(Aco, file=paste(fname, "score Aco_train_overlap.Robj", sep="/"))
  save(Aco_plus, file=paste(fname, "score Aco+_train_overlap.Robj", sep="/"))
  save(Bco, file=paste(fname, "score Bco_train_overlap.Robj", sep="/"))
  save(Bco_min, file=paste(fname, "score Bco-_train_overlapl.Robj", sep="/"))
  
  
}
