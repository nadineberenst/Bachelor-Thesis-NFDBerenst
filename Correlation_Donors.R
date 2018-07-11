####Load all necessities 

#install.packages("ggrepel")
#install.packages("pheatmap")
#install.packages("xlsx")

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


####Define donors
AB_cells <- pancreas@data.info[(pancreas@data.info$CT_TSNE == "Alpha" | pancreas@data.info$CT_TSNE == "Beta"),]

D1 <- rownames(AB_cells[AB_cells$Donor == "D1",])
D2 <- rownames(AB_cells[AB_cells$Donor == "D2",])
D3 <- rownames(AB_cells[AB_cells$Donor == "D3",])
D4 <- rownames(AB_cells[AB_cells$Donor == "D4",])
D6 <- rownames(AB_cells[AB_cells$Donor == "D6",])
D7 <- rownames(AB_cells[AB_cells$Donor == "D7",])
D8 <- rownames(AB_cells[AB_cells$Donor == "D8",])

####Calculate Pearson correlation
cordata<-cor(pancreas@data[,c(D1,D2,D3,D4,D6,D7,D8)])

annotation<-pancreas@data.info[c(D1,D2,D3,D4,D6,D7,D8),c("CT_TSNE","Donor")]


##Plot correlation in hetmap
ann_colors<-list("Donor"=c("D1"="lightcoral", "D2" = "midnightblue", "D3" = "#a6761d", "D4"="#FFFF00", "D6" = "darkmagenta", "D7"="deepskyblue", "D8"="bisque4"),
                 "CT_TSNE"=c("Alpha"="#7570b3","Beta"="#e7298a"))

plot.2.file("HM_Cor_allgenes_donors2")
pheatmap(cordata,kmeans_k = NA, breaks = NA, border_color = "grey60",cellwidth = NA,
         cellheight = NA, scale = "none", cluster_rows = TRUE,
         cluster_cols = TRUE, clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation", clustering_method = "ward.D2",  
         cutree_rows = NA, cutree_cols =  NA,
         treeheight_row = 0,
         treeheight_col = 40, legend = TRUE,
         legend_breaks = NA, legend_labels = NA, annotation_row = NA,
         annotation_col = annotation, annotation_colors = ann_colors,
         annotation_legend = TRUE, drop_levels = TRUE, show_rownames = F,
         show_colnames = F, main = NA, fontsize = 10, fontsize_row = 4,
         fontsize_col = 10, display_numbers = F, number_format = "%.2f",
         number_color = "grey30", fontsize_number = 0.8 * e,
         gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
         labels_col = NULL, filename = NA, width = NA, height = NA,     #ColSideColors=clab,
         silent = FALSE)
dev.off()
 