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

load("significantgenes_islets.Robj")
isletgenes <- genes
load("significantgenes_alphacells.Robj")
alphagenes <- genes
load("significantgenes_betacells.Robj")
betagenes <- genes

alphagenes <- (unique(sub("(\\_|\\.).+","", alphagenes)))
betagenes <- (unique(sub("(\\_|\\.).+","", betagenes)))
isletgenes <- (unique(sub("(\\_|\\.).+","", isletgenes)))

####Compare to set of bulk genes found in literature
#Segerstolpe et al
Alpha_Segerstolpe <- read.xlsx("Alphagenes_Segerstolpe.xlsx",1, as.data.frame=F)
Alpha_Segerstolpe <- as.character(Alpha_Segerstolpe)
Beta_Segerstolpe <- read.xlsx("Betagenes_Segerstolpe.xlsx", 1, as.data.frame=F)
Beta_Segerstolpe <- as.character(Beta_Segerstolpe)
Alpha_Xin <- read.xlsx("Alphagenes_Xin.xlsx", 1, as.data.frame=F)
Alpha_Xin <- as.character(Alpha_Xin)
Beta_Xin <- read.xlsx("Betagenes_Xin.xlsx", 1, as.data.frame=F)
Beta_Xin <- as.character(Beta_Xin)
Bulk_Xin <- read.xlsx("Bulkgenes_Xin.xlsx", 1, as.data.frame=F)
Bulk_Xin <- as.character(Bulk_Xin)
                            
Overlap_S_A <- calculate.overlap(list(Alpha_Segerstolpe, alphagenes))
Overlap_S_B <- calculate.overlap(list(Beta_Segerstolpe, betagenes))
Overlap_X_A <- calculate.overlap(list(Alpha_Xin, alphagenes))
Overlap_X_B <- calculate.overlap(list(Beta_Xin, betagenes))    
Overlap_X_Bulk <- calculate.overlap(list(Bulk_Xin, isletgenes))
Overlap_SX_A <- calculate.overlap(list(Alpha_Segerstolpe, Alpha_Xin))
