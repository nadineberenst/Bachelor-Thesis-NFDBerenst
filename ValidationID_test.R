####Load all necessities

setwd("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/ValidationID")

load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/pancreas.Robj")

source("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Required functions.R")

#install.packages("Metrics")
library(Metrics)
library(dplyr)
library(ggrepel)
library(robustbase)
library(Seurat)
library(base)


####Create test-loop
{
donors <- c("D1", "D2", "D3", "D4", "D6", "D7", "D8")
load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/score_all_donors.Robj")
score_all_7 <- score_all

healthy <- pancreas@data.info[(pancreas@data.info$Condition == "Healthy"),]
healthy <- grep("P546", rownames(healthy), invert=T, value=T) 
healthy <- pancreas@data.info[healthy,]  

for (i in 1:length(donors)){
  
  fname <- paste(donors[i], "test", sep="-")
  load(file=paste(fname, "DG", sep="/"))
  load(file=paste(fname, "score DG Alpha_train.Robj", sep="/"))
  load(file=paste(fname, "score DG Beta_train.Robj", sep="/"))
  load(file=paste(fname, "score Alpha mean_all.Robj", sep="/"))
  load(file=paste(fname, "score Beta mean_all.Robj", sep="/"))
  load(file=paste(fname, "score Aco_all.Robj", sep="/"))
  load(file=paste(fname, "score Aco+_all.Robj", sep="/"))
  load(file=paste(fname, "score Bco_all.Robj", sep="/"))
  load(file=paste(fname, "score Bco-_all.Robj", sep="/"))
  
  
  test_info <- healthy[grep(donors[i], healthy$Donor),]
  test_info <- na.omit(test_info)
  
  a <- rownames(test_info[(test_info$CT_TSNE == "Beta"),]); a2= "Beta"
  b <- rownames(test_info[(test_info$CT_TSNE == "Alpha"),]); b2 = "Alpha"
  
  
  test_data <- pancreas@data[,rownames(test_info)]
  pancdata<- test_data / rowMeans(test_data)
  
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
  Dtest <- score
  
  save(Dtest, file=paste(fname, "score_all_test.Robj", sep="/"))
  save(pancdata, file=paste(fname, "pancdata_test.Robj", sep="/"))
 
All7 <- score_all_7[names(Dtest)]
test_score <- cbind(All7, Dtest) 


save(test_score, file=paste(fname, "scores_7andtest.Robj", sep="/"))


All7M <- data.frame(All7)
DtestM <- data.frame(Dtest)
All7M$Condition <- "7Donors"
DtestM$Condition <- "Test donor"
names(DtestM) <- c("score_all", "Condition")
names(All7M) <- names(DtestM)
hist <- rbind(DtestM, All7M)
colors = c("7Donors" = "gold", "Test donor" = "blue4")

RMSE <- rmse(All7, Dtest)
RMSE  

plot.2.file(name=paste(fname, "histID", sep="/"))
p=ggplot(hist, aes(x=score_all, fill=Condition))+ geom_histogram(position="dodge")+
  scale_fill_manual(values=colors)+
  scale_x_continuous(limits = c(-0.7, 1.5))+
  theme_minimal(base_size=16) +
  ggtitle(paste(paste("IDscore","Validation",sep="-"), RMSE, sep=":"))+ xlab("ID-score") +ylab("Count")+ 
  geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
             linetype = "dashed", alpha=0.6)+ 
  geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
             linetype ="solid", alpha=0.6) 
print(p)
dev.off()



plot.2.file(name=paste(fname, "densID", sep="/"))
p=ggplot(hist, aes(x=score_all, color=Condition, y=..count..))+ geom_density()+
  scale_color_manual(values=colors)+
  scale_x_continuous(limits = c(-0.7, 1.5))+
  theme_minimal(base_size=16) +
  ggtitle(paste(paste("IDscore","Validation",sep="-"), RMSE, sep=":"))+ xlab("ID-score")+ 
  geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
             linetype = "dashed", alpha=0.6)+ 
  geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
             linetype ="solid", alpha=0.6) 
print(p)
dev.off()




save(RMSE, file=paste(fname, "rmse", sep="/"))
}

####Check RMSE
load("D1-test/rmse")
D1_rmse <- RMSE
load("D2-test/rmse")
D2_rmse <- RMSE
load("D3-test/rmse")
D3_rmse <- RMSE
load("D4-test/rmse")
D4_rmse <- RMSE
load("D6-test/rmse")
D6_rmse <- RMSE
load("D7-test/rmse")
D7_rmse <- RMSE
load("D8-test/rmse")
D8_rmse <- RMSE
}

####Calculate for overlapping DG's

donors <- c("D1", "D2", "D3", "D4", "D6", "D7", "D8")


healthy <- pancreas@data.info[(pancreas@data.info$Condition == "Healthy"),]
healthy <- grep("P546", rownames(healthy), invert=T, value=T) 
healthy <- pancreas@data.info[healthy,]  

for (i in 1:length(donors)){
  
  fname <- paste(donors[1], "test", sep="-")
  
  load(file=paste(fname, "DG_overlap", sep="/"))
  load(file=paste(fname, "score DG Alpha_train_overlap.Robj", sep="/"))
  load(file=paste(fname, "score DG Beta_train_overlap.Robj", sep="/"))
  load(file=paste(fname, "score Alpha mean_train_overlap.Robj", sep="/"))
  load(file=paste(fname, "score Beta mean_train_overlap.Robj", sep="/"))
  load(file=paste(fname, "score Aco_train_overlap.Robj", sep="/"))
  load(file=paste(fname, "score Aco+_train_overlap.Robj", sep="/"))
  load(file=paste(fname, "score Bco_train_overlap.Robj", sep="/"))
  load(file=paste(fname, "score Bco-_train_overlapl.Robj", sep="/"))
  
  genes <- overlap_genes
  
  test_info <- healthy[grep(donors[1], healthy$Donor),]
  test_info <- na.omit(test_info)
  
  a <- rownames(test_info[(test_info$CT_TSNE == "Beta"),]); a2= "Beta"
  b <- rownames(test_info[(test_info$CT_TSNE == "Alpha"),]); b2 = "Alpha"
  
  
  test_data <- pancreas@data[,rownames(test_info)]
  pancdata<- test_data / rowMeans(test_data)
  
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
  Dtest <- score
  
  save(Dtest, file=paste(fname, "score_all_test_overlap.Robj", sep="/"))
  save(pancdata, file=paste(fname, "pancdata_test_overlap.Robj", sep="/"))
  
  load("score_all7_overlap.Robj")
  All7 <- score_7overlap[names(Dtest)]
  test_score <- cbind(All7, Dtest) 
  
  
  save(test_score, file=paste(fname, "scores_7andtest_overlap.Robj", sep="/"))
  
  
  All7M <- data.frame(All7)
  DtestM <- data.frame(Dtest)
  All7M$Condition <- "7Donors"
  DtestM$Condition <- "Test donor"
  names(DtestM) <- c("score_all", "Condition")
  names(All7M) <- names(DtestM)
  hist <- rbind(DtestM, All7M)
  colors = c("7Donors" = "gold", "Test donor" = "blue4")
  
  plot.2.file(name=paste(fname, "histID", sep="/"))
  p=ggplot(hist, aes(x=score_all, fill=Condition))+ geom_histogram(position="dodge")+
    scale_fill_manual(values=colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("IDscore","Validation",sep="-"))+ xlab("ID-score") +ylab("Count")+ 
    geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
  
plot.2.file(name=paste(fname, "densID", sep="/"))
  p=ggplot(hist, aes(x=score_all, color=Condition, y=..count..))+ geom_density()+
    scale_color_manual(values=colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("IDscore","Validation",sep="-"))+ xlab("ID-score")+ 
    geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
  
  
  RMSE <- rmse(All7, Dtest)
 
  
  save(RMSE, file=paste(paste(fname, "rmse", sep="/"), "overlap", sep="_"))
       
title <- paste(fname, RMSE, sep=":")
       
plot.2.file(name=paste(fname, "histID_RMSE", sep="/"))
    p=ggplot(hist, aes(x=score_all, fill=Condition))+ geom_histogram(position="dodge")+
         scale_fill_manual(values=colors)+
         scale_x_continuous(limits = c(-0.7, 1.5))+
         theme_minimal(base_size=16) +
         ggtitle(title)+ xlab("ID-score") +ylab("Count")+ 
         geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
                    linetype = "dashed", alpha=0.6)+ 
         geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
                    linetype ="solid", alpha=0.6) 
       print(p)
       dev.off()
       
       
       
}

####Check RMSE
load("D1-test/rmse_overlap")
D1_rmse <- RMSE
load("D2-test/rmse_overlap")
D2_rmse <- RMSE
load("D3-test/rmse_overlap")
D3_rmse <- RMSE
load("D4-test/rmse_overlap")
D4_rmse <- RMSE
load("D6-test/rmse_overlap")
D6_rmse <- RMSE
load("D7-test/rmse_overlap")
D7_rmse <- RMSE
load("D8-test/rmse_overlap")
D8_rmse <- RMSE