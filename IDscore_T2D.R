#### Load necessities 
{
  setwd("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D")
  
  load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/pancreas.Robj")
  
  
  source("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Required functions.R")
  
  library(dplyr)
  library(ggrepel)
  library(robustbase)
  library(Seurat)
  library(base)
  
#Load signature genes, means and standard deviations calculated for 7 donors
  
  load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/score Aco_all.Robj")
  load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/score Bco_all.Robj")
  load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/score Alpha mean_all.Robj")
  load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/score Beta mean_all.Robj")
  load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/score DG Alpha_donors.Robj")
  load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/score DG Beta_donors.Robj")
  load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/score Aco+_all.Robj")
  load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/score Bco-_all.Robj")
  load("significantgenes.Robj")
    DG_7D <- genes
  
}

####Calculate scores for alpha and beta cells of diabetic subjects

#Pick alpha and beta cells
  patients <- pancreas@data.info[(pancreas@data.info$Condition == "T2D"),]
 
  Alpha_pat <- rownames(patients[patients$CT_TSNE =="Alpha",])
  Beta_pat <- rownames(patients[patients$CT_TSNE =="Beta",])
  
#Score cells 
  pancdata_pat <- pancreas@data[,rownames(patients)]
  pancdata_pat<-pancdata_pat/rowMeans(pancdata_pat)
  
  AB<-A-B
  AB<-na.omit(AB)
  C_pat<-pancdata_pat[DG_7D,c(Alpha_pat, Beta_pat)]
  CB_pat<-sweep(C_pat,1,B,"-");
  CB_pat<-na.omit(CB_pat)
  CB_pat<-sweep(C_pat,1,B,"-");
  t1<-sweep(CB_pat,1,AB,"*")
  t1<-colSums(t1)
  t2<-sqrt(sum(AB^2))
  score<-t1/(t2^2)
  score<-score[order(score, decreasing=F)]
  score_all <- score
  save(score_all, file="score_all_pat.Robj")
  
  
  mA<-score_all[score_all>Aco]
  mB<-score_all[score_all<Bco]
  Int<-score_all[score_all<Aco & score_all>Bco]
  length(Int)
 
####Annotate and plot
 scoreM<-data.frame(score_all)
  
 CellT<-pancreas@data.info[names(score_all),"CT_TSNE"]
 CellT <- as.data.frame(CellT, row.names=names(score_all))
 hist_pat<-cbind(scoreM,CellT)
 donors <- pancreas@data.info[names(score_all), "Donor"]
 donors <- as.data.frame(donors, row.names=names(score_all))
 hist_pat <- cbind(hist_pat, donors)
 names(hist_pat) <- c("score_all", "CellT", "Donor")
 hist_pat$Cells[(hist$CellT=="Alpha")] <- "Alpha T2D"
 hist_pat$Cells[(hist$CellT=="Alpha")] <- "Beta T2D"
  
#Plot ID score
  
  colors  = c("Alpha"= "#7570b3", "Beta"="#e7298a")
  
  plot.2.file(paste("T2D","dens_IDs_conditions_2stdev_7D",sep="_"))
  p=ggplot(hist_pat, aes(x=score_all, color=CellT y=..count..))+ geom_density()+
       scale_color_manual(values=colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("IDscore","T2D",sep="-"))+ xlab("ID-score") +ylab("Count")+ 
    geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
  
  plot.2.file(paste("T2D","hist_IDs_conditions_2stdev_7D",sep="_"))
  p=ggplot(hist_pat, aes(x=score_all, fill=CellT))+ geom_histogram()+
    scale_fill_manual(values=colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("IDscore","T2D",sep="-"))+ xlab("ID-score") +ylab("Count")+ 
    geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
  
  
  
  
  Colors= c("T2D1"= "#e6ab02", "T2D2"="#d95f02","T2D3"="#7570b3","T2D4"="#e7298a")
  
  plot.2.file(paste("T2D","hist_IDs_donors_2stdev_7D",sep="_"))
  p=ggplot(hist_pat, aes(x=score_all, fill=Donor))+ geom_histogram(alpha=0.6, position="stack")+    
    scale_fill_manual(values=Colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("7D_ID Donors","T2D",sep="-"))+ xlab("ID-score") +ylab("Count")+ 
    geom_vline(xintercept = c(Bco_min,Bco,Aco,Aco_plus), size = 1, colour =c("grey","grey", "grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
  
   
  plot.2.file(paste("T2D","dens_IDs_donors_2stdev_7D",sep="_"))
  p=ggplot(hist_pat, aes(x=score_all, color=Donor))+ geom_density(size=1)+    
    scale_color_manual(values=Colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("7D_ID Donors","T2D",sep="-"))+ xlab("ID-score") +
    geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
  }
#### Compare to ID-scores 7 healthy donors 
  
#Create variable with ID-scores of healthy donors 
  load("hist.Robj")
  hist$Cells[(hist$CellT=="Alpha")] <- "Alpha healthy"
  hist$Cells[(hist$CellT=="Beta")] <- "Beta healthy"
  hist_pat$Condition <- "T2D"
  hist$Condition <- "Healthy"
  data <- rbind(hist_pat, hist)
  
  
colors = c("Alpha healthy" = "grey", "Beta healthy"= "dark grey", "Alpha T2D" = "dark green", "Beta T2D" = "#66a61e")
plot.2.file(paste("T2D vs Healthy","hist_IDs_conditions_2stdev_7D",sep="_"))
    p=ggplot(data, aes(x=score_all, fill=Cells))+ geom_histogram(alpha=0.6)+    
    scale_fill_manual(values=colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("ID_7D","T2D",sep="-"))+ xlab("ID-score") +
    geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey", "grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
  
plot.2.file(paste("T2D vs Healthy","dens_IDs_conditions_2stdev_7D",sep="_"))
    p=ggplot(data, aes(x=score_all, color=Cells))+ geom_density(alpha=0.6)+    
    scale_color_manual(values=colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("ID_7D","T2D",sep="-"))+ xlab("ID-score") +
    geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey", "grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
  
Colors= c("T2D1"= "#e6ab02", "T2D2"="#d95f02","T2D3"="#7570b3","T2D4"="#e7298a", "D1"="lightcoral", "D2" = "midnightblue", "D3" = "#a6761d", "D4"="#FFFF00", "D6" = "darkmagenta", "D7"="deepskyblue", "D8"="bisque4")
Line = c("Healthy"= "dashed", "T2D" = "solid")


plot.2.file(paste("T2D vs Healthy","hist_IDs_donors_2stdev",sep="_"))
  p=ggplot(data, aes(x=score_all, fill=Donor))+ geom_histogram(alpha=0.6)+   
    scale_fill_manual(values=Colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("ID_7D","T2D",sep="-"))+ xlab("ID-score") +ylab("Count")+ 
    geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
 

plot.2.file(paste("T2D vs Healthy","dens_IDs_donors_2stdev",sep="_"))
  p=ggplot(data, aes(x=score_all, color=Donor, linetype=Condition, y=..count..))+ geom_density(size=1)+   
    scale_linetype_manual(values=Line)+
    scale_color_manual(values=Colors)+
    scale_x_continuous(limits = c(-0.7, 1.5))+
    scale_y_continuous(limits=c(0,2100))+
    theme_minimal(base_size=16) +
    ggtitle(paste("ID_7D","T2D",sep="-"))+ xlab("ID-score") +
    geom_vline(xintercept = c(Bco_min, Bco,Aco, Aco_plus), size = 1, colour =c("grey","grey","grey","grey"),
               linetype = "dashed", alpha=0.6)+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("grey","grey"),
               linetype ="solid", alpha=0.6) 
  print(p)
  dev.off()
  
  
####Save
 save(hist_pat, file="hist_pat.Robj")
 save(data, file="plot_donpat.Robj")
  

{ 
####Find T2D beta cells with a score higher than Bco 
 bc_pat <- rownames(hist_pat[hist_pat$CellT == "Beta",])
bc_score <- hist_pat[bc_pat,]
str_t2db <- rownames(bc_score[bc_score$score_all > Bco,]) 

save(str_t2db, file="strange T2D beta cells.Robj")}
