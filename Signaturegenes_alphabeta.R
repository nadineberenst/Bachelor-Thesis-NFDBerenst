####Load all necessities 

#install.packages("ggrepel")
#install.packages("pheatmap")

library(robustbase)
library(dplyr)
require(ggrepel)
require(pheatmap)

source("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Required functions.R")

setwd("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D")


load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/pancreas.Robj")

####Define alpha and beta cells 
healthy <- pancreas@data.info[(pancreas@data.info$Condition=="Healthy"),]

#Remove D5
Ialpha<-rownames(healthy[(healthy$CT_TSNE =="Alpha"),])
Ialpha <- grep("P546", Ialpha, invert=T, value=T)
Ibeta<-rownames(healthy[(healthy$CT_TSNE =="Beta"),])
Ibeta <- grep("P546", Ibeta, invert=T, value=T)

a<-Ibeta; a2<-"beta"; length(a)
b<-Ialpha; b2<-"alpha"; length(b)

A<-data.matrix(pancreas@data[,a])
B<-data.matrix(pancreas@data[,b])

####Perform Wilcoxon rank sum test
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
Wx2$max <- apply(Wx2[,c("Med_beta","Med_alpha")], 1, max) 

Wx2<-na.omit(Wx2)
Wx2<-Wx2[order(abs(Wx2$lfc2),decreasing = T),]

alphabeta_genes_ordered <- Wx2[Wx2$Padj<0.01 & abs(Wx2$lfc2)>1,]
alphabeta_top10 <- Wx2[1:10,]
genes<-rownames(Wx2[Wx2$Padj<0.01 & abs(Wx2$lfc2)>1,])

write.table(Wx2,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/Wx2.csv",quote=F,sep=";", row.names=T)
save(Wx2, file= "C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/Wx2")
write.table (alphabeta_top10, "Top 10 genes alpha vs beta.csv", row.names=T, quote=F, sep=";")
write.table(alphabeta_genes_ordered, "alphabeta_genes_ordered.csv", row.names=T, quote=F, sep=";")

#### Plot expression of DE genes
## Volcano plot
{
  VPd <- data.frame(row.names = row.names(Wx2),gene = row.names(Wx2),
                    pvalue = -log10(Wx2$Padj),
                    lfc = Wx2$lfc2)
  rownames(VPd)<-VPd$gene
  
  Wx2_up <- Wx2[order(Wx2$lfc2,decreasing = T),]
  Wx2_down <-Wx2[order(Wx2$lfc2,decreasing = F),]
  top_labelled<-VPd[c(rownames(Wx2_up[1:10,]), rownames(Wx2_down[1:10,])),]

  VPd <- VPd %>%
    mutate(color = ifelse(VPd$lfc > 1 & VPd$pvalue > -log10(0.01), 
                          yes = "P<.01, lfc>1",
                          no = ifelse(VPd$lfc < -1 & VPd$pvalue > -log10(0.01),
                                      yes = "P<.01, lfc<1", 
                                      no = "Rest")))
  
  
  ggtitle="Differentially expressed genes"
 # subtitle="All healthy donors"
  subtitle="7 healthy donors"
  
  
  simple <- ggplot(VPd, aes(x = lfc, y = pvalue)) + 
    # geom_point(size = 3, alpha = 0.7, na.rm = T) + # Make dots bigger
    geom_point(aes(col = factor(color)), size = 2, alpha = 0.7, na.rm = T) +
    scale_color_manual(values = c("P<.01, lfc>1" = "#e41a1c",
                                  "P<.01, lfc<1" = "#33a02c",
                                  "Rest" = "grey"))+ # change colors
    theme_bw(base_size = 16) + # change theme
    theme(legend.title = element_blank()) + # change theme
    ggtitle(label = ggtitle, subtitle = subtitle)+   # Add a title
    xlab(expression(log[2]("Beta" / "Alpha"))) + # x-axis label
    ylab(expression(-log[10]("adjusted p-value"))) + # y-axis label
    # geom_vline(xintercept = c(-2,2), colour = "darkgrey") + # Add cutoffs
    geom_hline(yintercept = -log10(0.01), colour = "darkgrey") + # Add cutoffs
    # geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
    # scale_colour_gradient(low = "black", high = "black", guide = FALSE) # Color black
    # scale_x_continuous(limits = c(-4, 4)) # min/max of lfc
    geom_text_repel(top_labelled,
                    mapping = aes(label = gene),
                    size = 2.5,
                    fontface = 1,
                    color = 'black',
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.5, "lines"))
  
simple 


plot.2.file(paste("SG_volcano_average", "alldonors", sep="-"), type="pdf")
plot.2.file(paste("SG_volcano_average", "7donors", sep="-"), type="pdf")
  simple
  dev.off()
  
  
  

#Volcano plot: log2 difference in median expression

  
  VPd <- data.frame(gene = row.names(Wx2),
                    pvalue = -log10(Wx2$Padj),
                    lfc = Wx2$lfc1)
  rownames(VPd)<-VPd$gene
  
 
  
  VPd <- VPd %>%
    mutate(color = ifelse(VPd$lfc > 1 & VPd$pvalue > -log10(0.01), 
                          yes = "P<.01, lfc>1",
                          no = ifelse(VPd$lfc < -1 & VPd$pvalue > -log10(0.01),
                                      yes = "P<.01, lfc<1", 
                                      no = "Rest")))
  
  
  simple <- ggplot(VPd, aes(x = lfc, y = pvalue)) + 
    # geom_point(size = 3, alpha = 0.7, na.rm = T) + # Make dots bigger
    geom_point(aes(col = factor(color)), size = 2, alpha = 0.7, na.rm = T) +
    scale_color_manual(values = c("P<.01, lfc>1" = "#e41a1c",
                                  "P<.01, lfc<1" = "#33a02c",
                                  "Rest" = "grey"))+ # change colors
    theme_bw(base_size = 16) + # change theme
    theme(legend.title = element_blank()) + # change theme
    ggtitle(label = ggtitle, subtitle = subtitle)+   # Add a title
    xlab(expression(log[2]("Beta" / "Alpha"))) + # x-axis label
    ylab(expression(-log[10]("adjusted p-value"))) + # y-axis label
    # geom_vline(xintercept = c(-2,2), colour = "darkgrey") + # Add cutoffs
    geom_hline(yintercept = -log10(0.01), colour = "darkgrey") + # Add cutoffs
    # geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
    # scale_colour_gradient(low = "black", high = "black", guide = FALSE) # Color black
    # scale_x_continuous(limits = c(-4, 4)) # min/max of lfc
    geom_text_repel(top_labelled,
                    mapping = aes(label = gene),
                    size = 2.5,
                    fontface = 1,
                    color = 'black',
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.5, "lines"))
  
  simple
  
#plot.2.file(paste("SG_volcano_median", "alldonors", sep="-"), type="pdf")
plot.2.file(paste("SG_volcano_median", "7donors", sep="-"), type="pdf")
  simple
  dev.off()

####Find signature genes

genesab <- Wx2[genes,]

write.table(genesab,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/alphabetasiggenes.csv", quote=F, row.names=T ,sep=";")

## Plot heatmap

names(genesab) <- c("names", "p_val", "Padj", "lfc1", "lfc2", "Med_beta", "Med_alpha", "Av_beta", "Av_alpha", "max", "genes", v12="Signature")
rownames(genesab) <- genesab$names

Colors=list("CellT"=c("Alpha" = "#7570b3", "Beta" = "#e7298a"), "Donor" =c("D1"="lightcoral", "D2" = "midnightblue", "D3" = "#a6761d", "D4"="#FFFF00", "D6" = "darkmagenta", "D7"="deepskyblue", "D8"="bisque4"))

data<-pancreas@data[rownames(genesab),c(Ibeta,Ialpha)]
data<-log10(data)
annotation<-pancreas@data.info[c(Ibeta,Ialpha),"CT_TSNE"]
annotation <- as.data.frame(annotation)
rownames(annotation) <- c(Ibeta, Ialpha)
donors <- pancreas@data.info[c(Ibeta, Ialpha), "Donor"]
donors <- as.data.frame(donors, row.names = c(Ibeta, Ialpha))
annotation <- cbind(annotation, donors)
names(annotation) <- c("CellT", "Donor")

plot.2.file(paste("Heatmap_alphabeta","7 Donors",sep="-"))
  pheatmap(data,kmeans_k = NA, breaks = NA, border_color = NA,cellwidth = NA,
         cellheight = NA, scale = "none", cluster_rows = TRUE,
         cluster_cols = FALSE, clustering_distance_rows ="euclidean", #"euclidean",
         clustering_distance_cols = "correlation", clustering_method = "ward.D2", #clustering_callback=res$tree_col,
         cutree_rows = NA, #cutree_cols =  k,
         treeheight_row = 40,
         treeheight_col = 40, legend = TRUE,
         legend_breaks = NA, legend_labels = NA, annotation_row = NA,
         annotation_col = annotation, annotation_colors = Colors,
         annotation_legend = TRUE, drop_levels = FALSE, show_rownames = F,
         show_colnames = F, main = NA, fontsize = 10, fontsize_row = 4,
         fontsize_col = 10, display_numbers = F, number_format = "%.2f",
         number_color = "grey30", fontsize_number = 0.8 * e,
         gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
         labels_col = NULL, filename = NA, width = NA, height = NA,     #ColSideColors=clab,
         silent = FALSE)

  dev.off()


save(genes, file="C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/significantgenes.Robj")
write.csv(genes,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/alphabetagenes.csv",quote=F,row.names=F)
write.table(genesab,"C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/alphabetasiggenes.csv", quote=F, row.names=T ,sep=";")

  