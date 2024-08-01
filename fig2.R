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

#fig2d
#data.combined is the seuratObj is this study
celltype=c("Cd4_Tnaive","Cd4_Tact_Cxcr3","Cd4_Tact_Tnf","Cd4_Tm","Tprolif",
           "Treg","Treg_resting","Treg_CTL",
           "Th1","Th1_act_Cxcr6","Th17","Th17_Il22",
           "Cd8_Tm")

celltype_col=c("#9DB4CB","#59518A","#39398A","#6494B3","#60AFB3",
      "#427256","#79992F", "#D9A833",
      "#AF1E1F","#5F0C0C","#F476BE","#962120",
      "#989898")
names(celltype_col)= celltype
pdf("fig2d.pdf",5,5)
DimPlot(data.combined, reduction = 'umap', label=T,raster=FALSE,cols =celltype_col )
dev.off()

#fig2e
library(scRNAtoolVis)
degs <- FindAllMarkers(object = data.combined, 
                       only.pos = F,
                       logfc.threshold = 0.25) 

library(scRNAtoolVis)
jjVolcano(diffData = degs,
          tile.col = celltype_col,
          size  = 3.5,
          topGeneN=5,
          fontface = 'italic',
          cluster.order=celltype,
          log2FC.cutoff=0.25,
          polar = T) +
  ylim(-8,10)
ggsave("fig2e.pdf",width = 12,height = 12)

#fig2f
pdf("dot_dc.pdf",9,4)
DotPlot(data.combined, 
                features = genes)+
  scale_color_gradientn(colours = colorRampPalette(c("blue", "yellow", "red3"))(99))+
  RotatedAxis2()+scale_x_discrete("")+scale_y_discrete("")
dev.off()

#fig2g
#sigScores is the pathway scores from VISION analysis
library(ggtern)
library(scales)
Activated=sigScores["Activated t cell",]
Inhibited=sigScores["Inhibited t cell",]
Naive=sigScores["Naive T cell",]
group=data.combined$celltype
data=data.frame(Activated,Inhibited,Naive,group)
data$group=as.factor(data$group)
ggtern(data=data, aes(x=Naive,
                    y=Activated,
                    z=Inhibited,color=group))+geom_point(aes(size=0.5), alpha=0.5)+ scale_color_manual(values =celltype_col)+
  theme(tern.panel.background = element_rect(fill = "white"),
         tern.panel.grid.minor = element_line(color = "gray90"), 
        tern.axis.arrow.show = TRUE, 
        tern.axis.arrow.T = element_line(color ='#FB7D80', size = 0.1), 
        tern.axis.arrow.L = element_line(color = 'grey30', size = 0.1),
        tern.axis.arrow.R = element_line(color = '#6FC7CF', size = 0.1),
        tern.axis.arrow.text.L = element_text(color = 'black'),  
        tern.axis.arrow.text.T = element_text(color = 'black'),
        tern.axis.arrow.text.R = element_text(color = 'black'),
        tern.axis.arrow.sep = 0.1, 
        tern.panel.grid.major.T = element_line(color = 'gray92', linetype = 1, linewidth = 0.8), 
        tern.panel.grid.major.L = element_line(color = 'gray92', linetype = 1, linewidth = 0.8),
        tern.panel.grid.major.R = element_line(color = 'gray92', linetype = 1, linewidth = 0.8),
        tern.panel.grid.minor.T = element_line(color = 'gray94', linetype = 1, linewidth = 0.8), 
        tern.panel.grid.minor.L = element_line(color = 'gray94', linetype = 1, linewidth = 0.8),
        tern.panel.grid.minor.R = element_line(color = 'gray94', linetype = 1, linewidth = 0.8),
        tern.axis.title.L = element_text(color = 'grey30', size = 11),
        tern.axis.title.T = element_text(color = '#FB7D80', size = 11),
        tern.axis.title.R = element_text(color = '#6FC7CF', size = 11),
        tern.axis.text.L = element_text(size = 17,face = 'bold'),
        tern.axis.text.R = element_text(size = 17,face = 'bold'),
        tern.axis.text.T = element_text(size = 17,face = 'bold'),
        tern.axis.vshift = 0.04,
        tern.axis.line.T = element_line(linewidth = 0.8),
        tern.axis.line.R = element_line(linewidth = 0.8),
        tern.axis.line.L = element_line(linewidth = 0.8))
ggsave('triangle_directlabels2.pdf',width = 6,height = 6)

#fig2h
#scenic input
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(SCENIC)
cell.info <- data.combined@meta.data
dge <- as.matrix(data.combined@assays$RNA@counts)
cenicOptions <- initializeScenic(
  org="mgi",   dbDir="./cisTarget_databases", # RcisTarget databases location
   dbs="mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
   datasetTitle="SCENIC on Mouse Cell Atlas", # choose a name for your analysis
  nCores=1
)
genesKept <- geneFiltering(
  exprMat=dge, 
  scenicOptions=scenicOptions,
  minCountsPerGene = 1,
  minSamples = 20
)
dge <- dge[genesKept, ]
dim(dge)
if(!dir.exists("output")) {
  dir.create("output")
}
saveRDS(cell.info, "output/s1_cell.info.rds")
source("utils/add_cellAnnotation.R")
saveLoom <- function(exprMat, output){
  ## prepare cell cluster info
  cellInfo2 <- data.frame(
    row.names = rownames(cell.info),
    #cellType = sapply(strsplit(colnames(exprMat), split = "\\."), function(x) x[1])
    cellType =cell.info$fincell
  )
  loom <- build_loom(output, dgem=exprMat)
  loom <- add_cellAnnotation(loom, cellInfo2)
  close_loom(loom)
}
saveLoom(dge, "output/s1_avg20_rep1.loom")
loom <- build_loom("output/s1_exprMat.loom", dgem=dge)
loom <- add_cellAnnotation(loom, cell.info)
close_loom(loom)

#Run scenic using python
#plot in R
rasMat <- readRDS("./output/s5_avg20_rep1.rasMat.rds")
s3_avg20_rep1.regulons <- read.delim("./output/s3_avg20_rep1.regulons.txt", header=FALSE, row.names=1)
rasMat2=t(rasMat)
rasMat_sub=rasMat2[,colnames(rasMat2)%in%rownames(data.combined@meta.data)]
row.names(rasMat_sub)=row.names(s3_avg20_rep1.regulons) 
sigScores <-as.matrix(rasMat_sub)
Idents(data.combined) <- "celltype"
data.combined@assays$RNA@counts=sigScores
data.combined@assays$RNA@data=sigScores
data.combined@assays$RNA@scale.data=sigScores
data.combined.markers <- FindAllMarkers(data.combined, 
                                        only.pos = TRUE, 
                                        min.pct = 0.25,
                                        assay ="RNA" ,
                                        logfc.threshold = 0)
saveRDS(data.combined.markers, "degs.rds")
top10 <- data.combined.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)

RotatedAxis2 <- function (...) 
{
  rotated.theme <- theme(axis.text.x = element_text(angle = 90, 
                                                    hjust = 1), validate = TRUE, ...)
  return(rotated.theme)
}

pdf("fig2h.pdf",9,4)
DotPlot(data.combined, 
                features = unique(top10$gene))+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  RotatedAxis2()+scale_x_discrete("")+scale_y_discrete("")
dev.off()

#fig2i&j
library(monocle)
RA_matrix<-as(as.matrix(data.combined@assays$RNA@counts), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(RA_matrix),gene_short_name=rownames(RA_matrix))
rownames(feature_ann)<-rownames(RA_matrix)
RA_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<- data.combined@meta.data
rownames(sample_ann)<-colnames(RA_matrix)
RA_pd<-new("AnnotatedDataFrame", data =sample_ann)
monocle_cds<-newCellDataSet(RA_matrix,phenoData =RA_pd,featureData =RA_fd,expressionFamily=negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
disp_table <- dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.2)
monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
monocle_cds <- reduceDimension(monocle_cds, reduction_method = 'tSNE', 
                              )
monocle_cds <- clusterCells(monocle_cds)
monocle_cds <- reduceDimension(monocle_cds,  reduction_method = "DDRTree",
                               num_dim = 5, scaling=T,
                               auto_param_selection = T
)
monocle_cds <- orderCells(monocle_cds)


p1 <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 1.5)+
  scale_color_gradientn(colors  = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))
p2 <- plot_cell_trajectory(monocle_cds, color_by = "clone", cell_size = 0.5)
cowplot::plot_grid(p1,p2, ncol = 2)
ggsave(paste0('fig2i_1','.pdf'), width=8, height=5)

library(ggridges)
plotdf=pData(monocle_cds)
plotdf$fincell=factor(plotdf$fincell,levels = celltype[celltype%in%unique(plotdf$fincell)])
pdf("fig2i_2.pdf",width = 5,height = 2.5)
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

library(patchwork)
library(reshape2)
col= c("grey90","#F05662")
data_plotC <- table(monocle_cds$clone, monocle_cds$fincell) %>% melt()
colnames(data_plotC) <- c("Sample", "CellType","Number")
data_plotC$CellType=factor(data_plotC$CellType,levels = celltype[celltype%in%unique(data_plotC$CellType)])
data_plotC$Sample=factor(data_plotC$Sample,levels = c("no","clone"))
pC1<-ggplot(data = data_plotC, aes(x = CellType, y = Number, fill = Sample)) +
  geom_bar(stat = "identity", width=0.8,aes(group=Sample),position="stack")+
  scale_fill_manual(values=col) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Average number")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) 
pC2<-ggplot(data = data_plotC, aes(x = CellType, y = Number, fill = Sample)) +
  geom_bar(stat = "identity", width=0.8,aes(group=Sample),position="fill")+
  scale_fill_manual(values=col) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
   theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))      
pC <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect')
ggsave(paste0("fig2i_3.pdf"), plot = pC, width = 7.5, height = 5)
