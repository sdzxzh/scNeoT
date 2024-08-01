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

#fig8a
#refer to fig2d

#fig8b
#refer to fig2e

#fig8c
#refer to fig4c

#fig8d
#refer to fig2g

#fig8e
sigScores <- as.matrix(t(getSignatureScores(vis)))
data=as.data.frame(t(sigScores))
type=factor(data.combined_fib@meta.data$celltype,levels=celltype)
group=as.vector(data.combined_fib@meta.data$group
my_comparisons <- list(c("NS", "NeoT"), 
                       c("DC", "NeoT"), 
                       c("OCDC", "NeoT"),
                       c("T", "NeoT") )
path=rownames(sigScores)
for (i in path) {
  Signature_score=data[i]
  Signature_score=Signature_score[,1]
  pdata_melt=data.frame(id,group,type,Signature_score)
  c <- ggplot(pdata_melt,
              aes(x=type, y=Signature_score, 
                  fill = group )) + 
                 geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = sample_color) +
       theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")
  p <- c + stat_compare_means(label = "p.signif")
  ggsave(paste0("cell&group_",i,".pdf"),plot =p, width=12, height=5)
  c <- ggplot(pdata_melt,
              aes(x=group, y=Signature_score, 
                  fill = group)) + 
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = sample_color) +
        theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")
  p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
  ggsave(paste0("group_",i,".pdf"),plot =p, width=2.4, height=5)
  c <- ggplot(pdata_melt,
              aes(x=type, y=Signature_score, 
                  fill = type, 
                  #color = group
              )) + 
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = celltype_col) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")
  p <- c + stat_compare_means(label = "p.signif")
  ggsave(paste0("cell_",i,".pdf"), plot =p,width=8, height=5)
  for (t in celltype) { 
    tmp <- pdata_melt[pdata_melt$type==t,]   
     c <- ggplot(tmp,
                aes(x=group, y=Signature_score, 
                    fill = group, 
                )) + 
      geom_boxplot(notch = F, alpha = 0.95, 
                   outlier.shape = 16,
                   outlier.colour = "black", 
                   outlier.size = 0.65) +
      scale_fill_manual(values = sample_color) +
           theme_classic() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
            axis.text.y = element_text(angle = 90, size = 12),
            axis.title.y = element_text(angle = 90, size = 15)) +
      theme(legend.position = "top")
    p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
    ggsave(paste0("split_cell_",i,"_",t,".pdf"), plot =p,width=2.4, height=5)
    
  }
}


#figS15a,b
sigScores <- as.matrix(t(getSignatureScores(vis)))
data=as.data.frame(t(sigScores))
type=factor(data.combined_fib@meta.data$celltype,levels=celltype)
group=as.vector(data.combined_fib@meta.data$group
my_comparisons <- list(c("NS", "NeoT"), 
                       c("DC", "NeoT"), 
                       c("OCDC", "NeoT"),
                       c("T", "NeoT") )
path=rownames(sigScores)
for (i in path) {
  Signature_score=data[i]
  Signature_score=Signature_score[,1]
  pdata_melt=data.frame(id,group,type,Signature_score)
  c <- ggplot(pdata_melt,
              aes(x=type, y=Signature_score, 
                  fill = group )) + 
                 geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = sample_color) +
       theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")
  p <- c + stat_compare_means(label = "p.signif")
  ggsave(paste0("cell&group_",i,".pdf"),plot =p, width=12, height=5)
  c <- ggplot(pdata_melt,
              aes(x=group, y=Signature_score, 
                  fill = group)) + 
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = sample_color) +
        theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")
  p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
  ggsave(paste0("group_",i,".pdf"),plot =p, width=2.4, height=5)
  c <- ggplot(pdata_melt,
              aes(x=type, y=Signature_score, 
                  fill = type, 
                  #color = group
              )) + 
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = celltype_col) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")
  p <- c + stat_compare_means(label = "p.signif")
  ggsave(paste0("cell_",i,".pdf"), plot =p,width=8, height=5)
  for (t in celltype) { 
    tmp <- pdata_melt[pdata_melt$type==t,]   
     c <- ggplot(tmp,
                aes(x=group, y=Signature_score, 
                    fill = group, 
                )) + 
      geom_boxplot(notch = F, alpha = 0.95, 
                   outlier.shape = 16,
                   outlier.colour = "black", 
                   outlier.size = 0.65) +
      scale_fill_manual(values = sample_color) +
           theme_classic() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
            axis.text.y = element_text(angle = 90, size = 12),
            axis.title.y = element_text(angle = 90, size = 15)) +
      theme(legend.position = "top")
    p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
    ggsave(paste0("split_cell_",i,"_",t,".pdf"), plot =p,width=2.4, height=5)
    
  }
}

#figS15c
library("RColorBrewer")
RotatedAxis2 <- function (...) 
{
  rotated.theme <- theme(axis.text.x = element_text(angle = 90, 
                                                    hjust = 1), validate = TRUE, ...)
  return(rotated.theme)
}
Idents(data.combined) <- "celltype"
pdf("fig7i_1.pdf",width = 12,height = 4.5)
DotPlot(data.combined_fib, features = gene)+RotatedAxis2()+
  scale_x_discrete("")+scale_y_discrete("") +scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))
dev.off()
Idents(data.combined_fib) <- "group"
pdf("figS15C_2.pdf",width = 12,height = 3)
DotPlot(data.combined_fib, features = gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("") +scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))
dev.off()




























