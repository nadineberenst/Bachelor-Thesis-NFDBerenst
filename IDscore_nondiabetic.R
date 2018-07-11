#### Load necessities 
{
source("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Required functions.R")
  
library(dplyr)
library(ggrepel)
library(robustbase)
library(Seurat)
library(base)
  
load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/pancreas.Robj")

setwd("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D")


load("significantgenes.Robj")

}

####ID Score computation
DG_all<- genes

healthy <- pancreas@data.info[(pancreas@data.info$Condition=="Healthy"),]
healthy <- grep("P546", rownames(healthy), invert=T, value=T)
healthy <- pancreas@data.info[healthy,]
  
Alphacells_all <-rownames(healthy[(healthy$CT_TSNE =="Alpha"),])
Betacells_all <- rownames(healthy[(healthy$CT_TSNE =="Beta"),])

save(Alphacells_all, file="Alphacells_all")
save(Betacells_all, file="Betacells_all")

pancdata <- pancreas@data[,rownames(healthy)]
pancdata_all<-pancdata/rowMeans(pancdata)
A<-rowMeans(pancdata_all[DG_all,Alphacells_all])
B<-rowMeans(pancdata_all[DG_all,Betacells_all])
C<-pancdata_all[DG_all,c(Alphacells_all, Betacells_all)]
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


save(score_all, file="score_all_donors.Robj")
save(pancdata_all, file="pancdata_all_donors.Robj")
save(A, file="score DG Alpha_donors.Robj")
save(B, file="score DG Beta_donors.Robj")

Amean<-mean(score_all[Alphacells_all])
Asd<-sd(score_all[Alphacells_all])
Aco<-Amean-2*Asd
Aco_plus <- Amean+2*Asd

Bmean<-mean(score_all[Betacells_all])
Bsd<-sd(score_all[Betacells_all])
Bco<-Bmean+2*Bsd
Bco_min <- Bmean-2*Bsd

save(Amean, file="score Alpha mean_all.Robj")
save(Bmean, file="score Beta mean_all.Robj")
save(Aco, file="score Aco_all.Robj")
save(Aco_plus, file="score Aco+_all.Robj")
save(Bco, file="score Bco_all.Robj")
save(Bco_min, file="score Bco-_all.Robj")

mA<-score_all[score_all>Aco]
mB<-score_all[score_all<Bco]
Int<-score_all[score_all<Aco & score_all>Bco]
length(Int)


scoreM<-data.frame(score_all)
CellT<-pancreas@data.info[names(score_all),"CT_TSNE"]
CellT<- as.data.frame(CellT, row.names=names(score_all))
hist<-cbind(scoreM,CellT)
donors <- pancreas@data.info[names(score_all), "Donor"]
donors <- as.data.frame(donors, row.names=names(score_all))
hist <- cbind(hist, donors)
names(hist) <- c("score_all", "CellT", "Donor")
Alpha = expression(alpha)

## plot IDscore
{
  plot.2.file(paste("ID_7","hist_IDs_conditions_2stdev",sep="_"))
  ggplot(hist, aes(x=score_all, fill = CellT)) + 
     scale_fill_manual(values=c("Alpha" = "#7570b3", "Beta"="#e7298a"))+
    geom_histogram(binwidth = 0.05, colour = "#1F3552") +
    scale_x_continuous(limits = c(-0.7, 1.5))+
    theme_minimal(base_size=16) +
    ggtitle(paste("ID_all","Celltype",sep="-"))+ xlab("ID-score") +ylab("Count")+ 
    geom_vline(xintercept = c(Bco_min,Bco,Aco, Aco_plus), size = 1, colour =c("#e7298a", "#e7298a","#7570b3", "#7570b3"),
               linetype = "dashed")+ 
    geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("#e7298a","#7570b3"),
               linetype ="solid")#+
  # annotate("text", x = 0.2:1.2, y = 70:70, label = c(sigma,sigma))
  dev.off()
  


plot.2.file(paste("ID_7","dens_IDs_conditions_2stdev",sep="_"))
ggplot(hist, aes(x=score_all, color=CellT)) + geom_density()+
  scale_color_manual(values=c("Alpha" = "#7570b3", "Beta"="#e7298a"))+
  scale_x_continuous(limits = c(-0.7, 1.5))+
  theme_minimal(base_size=16) +
  ggtitle(paste("ID_7","Celltype",sep="-"))+ xlab("ID-score") + 
  geom_vline(xintercept = c(Bco_min,Bco,Aco, Aco_plus), size = 1, colour =c("#e7298a", "#e7298a","#7570b3", "#7570b3"),
             linetype = "dashed")+ 
  geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("#e7298a","#7570b3"),
             linetype ="solid")#+
# annotate("text", x = 0.2:1.2, y = 70:70, label = c(sigma,sigma))
dev.off()

##Plot ID score per donor

title <- paste("ID_7","Donors",sep="-")


plot.2.file(paste("ID_7","dens_IDs_donors_2stdev",sep="_"))
ggplot(hist, aes(score_all, color = Donor)) + geom_density(alpha = 0.5) +
  scale_color_manual(values=c("D1"="lightcoral", "D2" = "midnightblue", "D3" = "#a6761d", "D4"="#FFFF00", "D5"="#336600", "D6" = "darkmagenta", "D7"="deepskyblue", "D8"="bisque4"))+
  scale_x_continuous(limits = c(-0.7, 1.5))+
  theme_minimal(base_size=16) +
  ggtitle(title)+ xlab("ID-score") +ylab("Count")+ 
  geom_vline(xintercept = c(Bco,Aco), size = 1, colour =c("#e7298a","#7570b3"),
             linetype = "dashed")+ 
  geom_vline(xintercept = c(Bmean,Amean), size = 1, colour =c("#e7298a","#7570b3"),
             linetype ="solid")#+
# annotate("text", x = 0.2:1.2, y = 70:70, label = c(sigma,sigma))
dev.off()
}

####Save scores of cells and histogram 

mB_beta<-mB[which(names(mB) %in% Betacells_all)]
mA_alpha<-mA[which(names(mA) %in% Alphacells_all)]
IDs = list(mB_betaD=mB_beta,mA_alpha=mA_alpha)

save(IDs, file="IDs_categorized.Robj")


save(hist, file= "hist.Robj")
