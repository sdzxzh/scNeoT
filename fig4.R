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


#fig4a
#data.combined is the seuratObj is this study
celltype=c("Cd8_Tnaive","Cd8_Tem","Cd8_Teff","Tex",
           "Cd4_Tnaive","Cd4_Tm","Cd4_Treg","Cd4_Th2","Tprolif","Cd4_Th1","Cd4_ISG","Cd4_Tanergic", "Cd4_THsp",
           "NK","NK_Xcl1")
celltype_col=c("#79992F","#D9A833","#AF1E1F","#F476BE","#B96C27",
      "#427256","#60AFB3","#BBD4E0","#6494B3","#39398A","#73398D","#462672","grey60","#511818",
      "#C37E7E","#5F0C0C"
)
names(celltype_col)= celltype
pdf("fig4a.pdf",6.1,4.6)
DimPlot(data.combined, reduction = 'umap', label=T,raster=FALSE,cols =celltype_col ) 
dev.off()

#fig4b
pdf("Atp2b4_2.pdf",9,3)
FeaturePlot(data.combined, raster=F,
            cols =c("lightgrey", "#D24641"),features =c("gene_name")  )
dev.off()

#fig4c
data_plotC <- table(data.combined@meta.data$group, data.combined@meta.data$celltype) %>% melt()
colnames(data_plotC) <- c("Sample", "CellType","Number")
pC1 <- ggplot(data = data_plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="stack")+
  scale_fill_manual(values=celltype_col) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Average number")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) 
pC2 <- ggplot(data = data_plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=celltype_col) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  #scale_y_continuous(labels = percent)+  
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))     
pC <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect')
ggsave(paste0("fig4c.pdf"), plot = pC, width = 6.5, height = 5)

#fig4d

heatmap_gene <- c("Tcf7", "Lef1","Il7r",
                  "Gzma", "Gzmb","Gzmk", "Eomes","IFNG"
                  "Pdcd1","Runx2",
                  "Foxp3","Ctla4", "Gata3",
                  "Mki67", "Top2a", 
                  "Cxcr3","Cxcr6","Ifit1","Ifit3",
                  "Tnfsf8","Nr4a1",
                  "Dnajb1","Hspa1a",
                  "Ccl5","Xcl1"
)
data.combined@meta.data$celltype <- factor(data.combined@meta.data$fincell, levels = celltype)
Idents(data.combined)="celltype"
heatmap_AveE <- AverageExpression(data.combined, assays = "RNA", features = heatmap_gene,verbose = TRUE) %>% .$RNA
gene_num <- c(3,5,2,3,2,4,2,2,2)
gaps_row <- cumsum(gene_num)
cluster_num <- c(rep(1,length(celltype)))
gaps_col <- cumsum(cluster_num)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
annotation_row <- data.frame(row.names = rownames(heatmap_AveE),
                             `GeneType` = rep(as.factor(c(1:length(gene_num))),gene_num))
annotation_col <- data.frame(row.names = colnames(heatmap_AveE),
                             `CellType` = colnames(heatmap_AveE))
annotation_colors = list(`CellType` = celltype_col)
pdf("fig4d_1.pdf",8,5)
pheatmap(heatmap_AveE,cluster_cols = F,cluster_rows = F,show_colnames=F,show_rownames=T,
                border=F,#border_color = "white",
                color = c(colorRampPalette(colors = c("#2166ac","#f7fbff"))(length(bk)/2),
                          colorRampPalette(colors = c("#f7fbff","#b2182b"))(length(bk)/2)),
                breaks=bk,scale="row",legend_breaks=seq(-2,2,2),
                gaps_row = gaps_row,
                gaps_col = gaps_col,
                #annotation_row = annotation_row,
                annotation_col = annotation_col,
                annotation_colors = annotation_colors,
                annotation_names_row = F,annotation_names_col = T)
dev.off()

library(reshape2)
library(ggpubr)
pB2_df <- table(data.combined@meta.data$celltype,data.combined@meta.data$group) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
pB2_df$Cluster <- factor(pB2_df$Cluster,levels = celltype)
for (i in unique(pB2_df$Sample)) {
  pB2_df[pB2_df$Sample==i,]$Number=pB2_df[pB2_df$Sample==i,]$Number/sum(pB2_df[pB2_df$Sample==i,]$Number)
}
sample_color <- c("#A48AD3","#1CC5FE","#6FC7CF","#FBA27D","#FB7D80")

pB2 <- ggplot(data = pB2_df, aes(x = Cluster, y = Number, fill = Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  rotate_x_text()+
  #scale_y_continuous(position = "right",labels = percent)+
  scale_y_continuous(position = "right")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))
stat_raw=data.frame(data.combined$celltype,data.combined$time)
table1=table(stat_raw$data.combined.celltype,stat_raw$data.combined.time)
table2=table1[,-5]
num_neot=sum(table1[,5])
num_other= sum(table2)  
ratio_neo=table1[,5]/num_neot
ratio_other=table2/num_other
fin_ratio=ratio_neo/ratio_other
AveExpression=data.frame(fin_ratio)
AveExpression$Gene <- "fin_ratio"
AveExpression$Cluster=rownames(AveExpression)
Ave_df=AveExpression
colnames(Ave_df) <- c("Expression", "Gene", "Cluster")
Ave_df$Group <- paste(Ave_df$Cluster,Ave_df$Gene,sep = "_")
Ave_df$Expression[which(Ave_df$Expression>3)] <- 3 Ave_df$Cluster=factor(Ave_df$Cluster,levels = celltype)
pB1 <- ggplot(Ave_df,aes(x=Gene, y=Expression))+
  geom_hline(yintercept = seq(0, 2, 0.5),linetype = 2, color = "lightgray",size=1)+
  geom_line()+
  geom_segment(aes(x=Gene,xend=Gene,y=0,yend=Expression),color="lightgray",size = 1.5)+
  geom_point(size=4,aes(color=Gene))+
  scale_color_manual(values=c("#FB7D80","#AF1E1F")) +
 # scale_color_continuous(values=c("#00AFBB", "#E7B800", "#FC4E07", "#41ab5d"))+
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="ratio")
pB1 <- facet(pB1, facet.by = "Cluster",ncol = length(unique(Ave_df$Cluster)),panel.labs.font = list(size = 12),panel.labs.background = list(fill = "#a6cee3"))
pB1 <- pB1 + scale_y_continuous(position = "right")+ 
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_blank())+ 
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "right",
        panel.border = element_blank(),
        axis.ticks.x = element_line(color =  NA))
require(ggplotify)
pB_1_2 <- pB1 + pB2 + plot_layout(ncol = 1, heights = c(1.2, 2))
ggsave(paste0("fig4d_2.pdf"), plot = pB_1_2, width = 8, height =5)

#fig4ef
#The construction method of the monocle object  is similar to fig2i
p1 <- plot_cell_trajectory(monocle_cds, color_by = "celltype", cell_size = 1.5)+
  scale_color_gradientn(colors  = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))
p2 <- plot_cell_trajectory(monocle_cds, color_by = "group", cell_size = 1.5)
cowplot::plot_grid(p1,p2, ncol = 2)
ggsave("fig4ef.pdf", width=8, height=5)

#fig4g
pseudotime_de <- differentialGeneTest(monocle_cds,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
pdf("fig4g_1.pdf",width = 4,height = 5)
plot_pseudotime_heatmap(monocle_cds[ pseudotime_de$gene_id[1:100], ], num_clusters = 3, cores=10, show_rownames=TRUE,
                        hmcols =colorRampPalette(rev(brewer.pal(9, "PRGn")))(100))
dev.off()
library(ggridges)
plotdf=pData(monocle_cds)
plotdf$celltype=factor(plotdf$celltype,levels = celltype)
pdf("fig4g_2.pdf",width = 5,height = 2.5)
ggplot(plotdf, aes(x=Pseudotime,y=fincell,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  #geom_vline(xintercept = c(5,10),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )

dev.off()

#fig4h
gene=c("Sell","Tcf7","Gzmk","Gzmb","Pdcd1","Ctla4")
plot_genes_in_pseudotime(monocle_cds[gene,],color_by = "fincell",cell_size = 2)+
  scale_colour_manual(values = celltype_col)
ggsave("fig4h.pdf", width=8, height=5)

#fig4i
library(ggtern)
library(scales)
library(ggplot2)
col=c("#A48AD3","#FB7D80","#6FC7CF","grey85","grey30")
names(col)=c("a_NS","e_NeOT","b_other_vacc","nosig","Mix")
ggtern(data=df, aes(x=a_NS,
                    y=b_other_vacc,
                    z=e_NeOT,color=group))+geom_point(aes(size=size), alpha=0.5)+ scale_color_manual(values =col)+
 # stat_density_tern(n = 100,h=1.2, alpha=0.7)+
  theme(tern.panel.background = element_rect(fill = "white"), 
        tern.panel.grid.minor = element_line(color = "gray90"), 
        tern.axis.arrow.show = TRUE, 
        tern.axis.arrow.T = element_line(color ='#A48AD3', size = 3), 
        tern.axis.arrow.L = element_line(color = '#6FC7CF', size = 3),
        tern.axis.arrow.R = element_line(color = '#FB7D80', size = 3),
        tern.axis.arrow.text.L = element_text(color = 'black'),  
        tern.axis.arrow.text.T = element_text(color = 'black'),
        tern.axis.arrow.text.R = element_text(color = 'black'),
        tern.axis.arrow.sep = 0.1, 
        tern.panel.grid.major.T = element_line(color = 'gray92', linetype = 1, size = 0.8), 
        tern.panel.grid.major.L = element_line(color = 'gray92', linetype = 1, size = 0.8),
        tern.panel.grid.major.R = element_line(color = 'gray92', linetype = 1, size = 0.8),
        tern.panel.grid.minor.T = element_line(color = 'gray94', linetype = 1, size = 0.8), 
        tern.panel.grid.minor.L = element_line(color = 'gray94', linetype = 1, size = 0.8),
        tern.panel.grid.minor.R = element_line(color = 'gray94', linetype = 1, size = 0.8),
        tern.axis.title.L = element_text(color = '#A48AD3', size = 11),
        tern.axis.title.T = element_text(color = '#6FC7CF', size = 11),
        tern.axis.title.R = element_text(color = '#FB7D80', size = 11),
        tern.axis.text.L = element_text(size = 17,face = 'bold'),
        tern.axis.text.R = element_text(size = 17,face = 'bold'),
        tern.axis.text.T = element_text(size = 17,face = 'bold'),
        tern.axis.vshift = 0.04,
        tern.axis.line.T = element_line(size = 0.8),
        tern.axis.line.R = element_line(size = 0.8),
        tern.axis.line.L = element_line(size = 0.8))
ggsave('fig4i.pdf')




















