####Load all necessities 

#install.packages("ggrepel")
#install.packages("pheatmap")

library(robustbase)
library(dplyr)
require(ggrepel)
require(pheatmap)
require(VennDiagram)
library(base)
require(gplots)

source("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Required functions.R")

#setwd("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset")
setwd("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D")


load("C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/pancreas.Robj")

####Define Alpha cells

patients <- pancreas@data.info[(pancreas@data.info$Condition == "T2D"),]


healthy <- pancreas@data.info[(pancreas@data.info$Condition == "Healthy"),]
healthy <- grep("P546", rownames(healthy), invert=T, value=T) 
healthy <- pancreas@data.info[healthy,]  


alpha_pat <- rownames(patients[(patients$CT_TSNE=="Alpha"),])
alpha_7D <- rownames(healthy[(healthy$CT_TSNE == "Alpha"),])


a<-alpha_pat; a2<-"alpha patients"; length(a)
b<-alpha_7D; b2<-"alpha healthy"; length(b)


A<-data.matrix(pancreas@data[,a]); dim(A)
B<-data.matrix(pancreas@data[,b]); dim(B)


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


Wx2<-na.omit(Wx2)
Wx2<-Wx2[order(abs(Wx2$lfc2),decreasing = T),]
alpha_genes_ordered <- Wx2[Wx2$Padj<0.01 & abs(Wx2$lfc2)>1,]

alpha_top10 <- Wx2[Wx2$Padj<0.01 & abs(Wx2$lfc2)>1,]
alpha_top10 <- alpha_top10[order(abs(alpha_top10$lfc2),decreasing = T),]
alpha_top10 <- alpha_top10[1:10,]

genes<-rownames(Wx2[Wx2$Padj<0.01 & abs(Wx2$lfc2)>1,])


write.table(Wx2,"alphacells_Wx2.csv",quote=F, row.names=T, sep=";")
save(Wx2, file= "alphacells_Wx2")
write.table(alpha_top10, "Top 10 alpha genes.csv", quote=F, row.names=T, sep=";")
write.table(alpha_genes_ordered, "Alphagenes_ordered.csv", quote=F, row.names=T, sep=";" )

#### Plot differentially expressed genes
## Volcano Plot

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
  subtitle="Alpha cells T2D vs healthy"
  
  simple <- ggplot(VPd, aes(x = lfc, y = pvalue)) + 
    # geom_point(size = 3, alpha = 0.7, na.rm = T) + # Make dots bigger
    geom_point(aes(col = factor(color)), size = 2, alpha = 0.7, na.rm = T) +
    scale_color_manual(values = c("P<.01, lfc>1" = "#e41a1c",
                                  "P<.01, lfc<1" = "#33a02c",
                                  "Rest" = "grey"))+ # change colors
    theme_bw(base_size = 16) + # change theme
    theme(legend.title = element_blank()) + # change theme
    ggtitle(label = ggtitle, subtitle = subtitle)+   # Add a title
    xlab(expression(log[2]("T2D" / "Healthy"))) + # x-axis label
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
  
  
  plot.2.file(paste("SG_volcano_average", "alpha cells", sep="-"), type="pdf")
  simple
  dev.off()
  
  
  #Volcano plot: log2 median difference
  
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
    xlab(expression(log[2]("T2D" / "Healthy"))) + # x-axis label
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
  
  plot.2.file(paste("SG_volcano_median", "alpha cells", sep="-"), type="pdf")
  simple
  dev.off()
  

## Create heatmap
Colors=list("Condition"=c("Healthy" = "grey", "T2D" = "#66a61e"), "Donor" =c("T2D1"= "#e6ab02", "T2D2"="#d95f02","T2D3"="#7570b3","T2D4"="#e7298a", "D1"="lightcoral", "D2" = "midnightblue", "D3" = "#a6761d", "D4"="#FFFF00", "D6" = "darkmagenta", "D7"="deepskyblue", "D8"="bisque4"))
  
data<-cbind(pancreas@data[genes,c(alpha_7D, alpha_pat)])
rownames(data) <- genes
data<-log10(data)
annotation<-c(rep("Healthy",length(alpha_8D)), rep("T2D", length(alpha_pat)))
annotation <- as.data.frame(annotation)
rownames(annotation) <- c(alpha_7D, alpha_pat)
names(annotation) <- "Condition"
Donor <- pancreas@data.info[c(alpha_7D, alpha_pat), "Donor"]
  
annotation<- cbind(annotation, Donor)
  
plot.2.file(paste("Heatmap", "Alpha T2D vs Healthy",sep="-"))
  pheatmap(data,kmeans_k = NA, breaks = NA, border_color = NA,cellwidth = NA,
           cellheight = NA, scale = "none", cluster_rows = T,
           cluster_cols = FALSE, clustering_distance_rows ="euclidean", #"euclidean",
           clustering_distance_cols = "correlation", clustering_method = "ward.D2", #clustering_callback=res$tree_col,
           cutree_rows = NA, #cutree_cols =  k,
           treeheight_row = 40,
           treeheight_col = 40, legend = TRUE,
           legend_breaks = NA, legend_labels = NA, annotation_row = NA,
           annotation_col = annotation, annotation_colors = Colors,
           annotation_legend = TRUE, drop_levels = FALSE, show_rownames = T,
           show_colnames = F, main = NA, fontsize = 10, fontsize_row = 6,
           fontsize_col = 10, display_numbers = F, number_format = "%.2f",
           number_color = "grey30", fontsize_number = 0.8 * e,
           gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
           labels_col = NULL, filename = NA, width = NA, height = NA,     #ColSideColors=clab,
           silent = FALSE)
  dev.off()
  

  
####Save signature genes    
save(genes, file="significantgenes_alphacells.Robj")
  
alphagenes <- Wx2[genes,]  
rownames(alphagenes) <- genes
short_alphagenes <- (unique(sub("(\\_|\\.).+","", rownames(alphagenes))))
short_alphagenes <- as.data.frame (short_alphagenes, row.names=rownames(alphagenes) )
alphagenes<- cbind(alphagenes, short_alphagenes)

write.table(alphagenes, file="C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/alphagenes.csv",  quote=F,row.names=T, sep=";")
  
sigalpha <- (unique(sub("(\\_|\\.).+","", genes)))

write.csv(sigalpha, file="C:/Users/Nadine/Documents/BEP/R directory/Seurat/Nathalie's code/Validation/Downsampling/T2 diabetes/One dataset/7D/alphagenes_short2.csv",  quote=F,row.names=F)
