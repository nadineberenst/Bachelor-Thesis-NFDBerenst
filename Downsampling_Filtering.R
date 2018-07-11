####Load all necessities

source("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Required functions.R")
library(devtools)
library(Seurat)
library(ggplot2)

setwd("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D")

load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/pancreas.Robj")  

####Remove donor 5 (P546)
cells_7D <- grep("P546", rownames(pancreas@data.info), invert=T, value=T)
tSNE_7D <- pancreas@tsne.rot[cells_7D,]
info_7D <- pancreas@data.info[cells_7D,]
data_7D <- pancreas@data[,cells_7D]
rawdata_7D <- pancreas@raw.data[,cells_7D]

####Check quality of dataset 
#####Create boxplots to check quality
healthy_7D <- info_7D[info_7D$Condition == "Healthy",]
patients_7D <- info_7D[info_7D$Condition == "T2D",]
quality_healthy <- healthy_7D[, c("nUMI", "nGene")]
quality_pat <- patients_7D[, c("nUMI", "nGene")]

quality_pat$Condition <- "T2D"
quality_healthy$Condition <- "Healthy"

data_all <- rbind(quality_pat, quality_healthy)
bp_UMI <- data_all[,c("nUMI", "Condition")]
bp_nGene <- data_all[,c("nGene", "Condition")]

Colors = c("Healthy" = "grey", "T2D" = "#66a61e")
Title = "Islet cells T2D vs healthy"

plot.2.file("bp_nUMI_data", type= "pdf")
ggplot(bp_UMI, aes(x=Condition, y=nUMI))+
  geom_boxplot(aes(fill=Condition))+
  scale_fill_manual(values=Colors)+
  ggtitle(Title)
dev.off()

plot.2.file("bp_nGene_data", type= "pdf")
ggplot(bp_nGene, aes(x=Condition, y=nGene))+
  geom_boxplot(aes(fill=Condition))+
  scale_fill_manual(values=Colors)+
  ggtitle(Title)
dev.off()


####tSNE visualization
sample<-109  #find for which sample the tSNE looks best

data.plot <- cbind(tSNE_7D, info_7D)

##Visualization donor dependence 
variable= "Donor"
Colors= c("T2D1"= "#e6ab02", "T2D2"="#d95f02","T2D3"="#7570b3","T2D4"="#e7298a", "D1"="lightcoral", "D2" = "midnightblue", "D3" = "#a6761d", "D4"="#FFFF00", "D6" = "darkmagenta", "D7"="deepskyblue", "D8"="bisque4")
position=c(0.8,0.8) # position legend
base=16 ## txt size graph
title=("tSNE")
subtitle=paste(variable, "distribution",sep=" ")

tsne<-ggplot(data.plot, aes(x=tSNE_1,y=tSNE_2))+
  geom_point(aes(colour=factor(Donor)),size=1.5)+   #can be changed to other variables
  #scale_color_manual(values= Colors) #to set colors manually (as defined in line 100)
  theme_bw(base_size = base)+ #adds grids to plot and scales it
  ggtitle(label = title, subtitle = subtitle) + labs(colour =paste(variable))+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust=0.5))+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position =position) #remove grid and axes, move legend to left side of the graph
tsne


plot.2.file(paste("tSNE_7donors_4pat", variable, sep=""))
tsne
dev.off()


##Plotting expression marker genes

Markergenes <- c("GCG__chr2","TTR__chr18","INS__chr11","SST__chr3","PPY__chr17","KRT19__chr17","PRSS1__chr7","COL1A1__chr17","CD24__chrY") 

dir.create("Markergenes_exp_norm", showWarnings = FALSE)
for (g in 1:length(Markergenes)){
  tryCatch({
      fname<-paste(Markergenes[g],"_tSNE overlay",sep="")
    plot.2.file(paste("Markergenes_exp_norm",fname,sep="/"),type="pdf")
    p=ggplot(tSNE_7D, aes(tSNE_1, tSNE_2),ax=FALSE, cex.lab=1.5,cex.main=1.5) +
      geom_point(aes(colour=t(data_7D[Markergenes[g], rownames(tSNE_7D)]))) +
      scale_colour_gradientn(colours=c(rev(terrain.colors(100))))+
      theme_bw(base_size = base)+ 
      ggtitle(paste(chop_chr(Markergenes[g])))+ labs(colour = "Transcripts") +
      theme(plot.title = element_text(hjust = 0.5))+ # , plot.subtitle = element_text(hjust=0.5)
      theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
            panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(legend.position =position,legend.title=element_text(size=base-4))
    print(p)
    dev.off()
  }, error=function(e){})
} 


#Raw expression of markergenes
dir.create("Markergenes_exp_raw", showWarnings = F)
for (g in 1:length(Markergenes)){
  tryCatch({
    fname<-paste(Markergenes[g],"_tSNE overlay",sep="")
    plot.2.file(paste("Markergenes_exp_raw",fname,sep="/"),type="pdf")
    p=ggplot(tSNE_7D, aes(tSNE_1, tSNE_2),ax=FALSE, cex.lab=1.5,cex.main=1.5) +
      geom_point(aes(colour=t(rawdata_7D[Markergenes[g], rownames(tSNE_7D)]))) +
      scale_colour_gradientn(colours=c(rev(terrain.colors(1000))))+
      theme_bw(base_size = base)+ 
      ggtitle(paste(chop_chr(Markergenes[g])))+ labs(colour = "Transcripts") +
      theme(plot.title = element_text(hjust = 0.5))+ # , plot.subtitle = element_text(hjust=0.5)
      theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
            panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(legend.position =position,legend.title=element_text(size=base-4))
    print(p)
    dev.off()
  }, error=function(e){})
} 


####Assign cell types to clusters

base=16 ## txt size graph
position = c(0.10, 0.70)
title=("tSNE map")
subtitle="Cell types"
variable="Cell types"

tsne<-ggplot(data.plot, aes(x=tSNE_1,y=tSNE_2))+
  geom_point(aes(colour=factor(CT_TSNE)),size=1.5)+   #can be changed to other variables
  #scale_color_manual(values= Colors) #to set colors manually (as defined in line 100)
  theme_bw(base_size = base)+ #adds grids to plot and scales it
  ggtitle(label = title, subtitle = subtitle) + labs(colour =paste(variable))+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust=0.5))+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position =position) #remove grid and axes, move legend to left side of the graph
tsne

plot.2.file(paste("tSNE_7donors_4pat", variable, sep=""))
tsne
dev.off()
