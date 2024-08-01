library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

#figS12b,c
#refer to fig4a,b

#figS12d
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
infercnv_obj = readRDS("./run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$cCON
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- test_loc$test
anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)
gn <- rownames(expr)
geneFile <- read.table("./geneFile.txt", header=FALSE, row.names=NULL)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
expr=expr[,colnames(expr)[test_loc]]#or normal_loc
set.seed(12345)
kmeans.result <- kmeans(t(expr), 7)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") 
kmeans_df_s=arrange(kmeans_df,kmeans_class) 
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) 
head(kmeans_df_s)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:19,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7] 
names(color_v)=as.character(1:7)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("a_NS"="#A48AD3",
                                                             "b_DC" = "#1CC5FE",
                                                             "c_OCDC" = "#6FC7CF",
                                                             "d_T" = "#FBA27D",
                                                             "e_NeOT" = "#FB7D80"
                                                             ),
                                                     kmeans_class=color_v))
pdf("figS12d.pdf",width = 15,height = 10)
ht = Heatmap(t(expr)[rownames(kmeans_df_s),], 
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), 
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:19,sep = "")), 
             column_gap = unit(2, "mm"),
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             top_annotation = top_anno,left_annotation = left_anno, 
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

#figS12e,f
#refer to fig4c,d

#figS12g
library(data.table)
library(pbapply)
library(plyr)
library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)
rasMat <- fread("output/s3_avg20_rep1.AUCell.txt", sep = "\t", header = T, data.table = F) # data.frame
rownames(rasMat) <- rasMat$V1
colnames(rasMat) <- sub("(+)", "", colnames(rasMat), fixed = T)
rasMat <- rasMat[, -1]
saveRDS(rasMat, "output/s5_avg20_rep1.rasMat.rds")
s3_avg20_rep1.regulons <- read.delim("./output/s3_avg20_rep1.regulons.txt", header=FALSE, row.names=1)
colnames(rasMat)=rownames(s3_avg20_rep1.regulons)
cell.info <- readRDS("output/s4_cell.info.rds")
cell.info$Tissue=cell.info$celltpye
cell.types <- names(table(cell.info$CellType))
ctMat <- lapply(cell.types, function(i) {
  as.numeric(cell.info$CellType == i)
})
ctMat <- do.call(cbind, ctMat)
colnames(ctMat) <- cell.types
rownames(ctMat) <- rownames(cell.info)
rssMat <- pblapply(colnames(rasMat), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
  })
})
rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(rasMat)
colnames(rssMat) <- colnames(ctMat)
saveRDS(rssMat, "output/s5_avg20_rep1.rssMat.rds")
rssMat <- readRDS("output/s5_avg20_rep1.rssMat.rds")
binMat <- read.table("output/s3_avg20_rep1.binary_mtx.txt", sep = "\t", header = T, row.names = 1, check.names = FALSE)
colnames(binMat) <- sub("(+)", "", colnames(binMat), fixed = T)
source("utils/plotRegulonRank.R")
PlotRegulonRank(rssMat, "the cell type we need", topn=5)







