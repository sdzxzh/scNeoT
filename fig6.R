library(Seurat)

#fig6b-e
#refer to fig5f-j

#fig6f
pdf("fig6f.pdf",20,3)
FeaturePlot(data.combined, 
                      cols=c('grey',"#CC3333"),
                        features =c("NeoT"))

dev.off()

#fig6g
df=table(seurat_neot$celltype)
df=df/cell_count
colnames(df)<-c("celltype","freq")
ggplot(df, aes(x = 1, y=freq, fill = celltype)) + 
  geom_bar(width = 1, stat = "identity")+
  xlim(1,2.5)+
  geom_text(aes(label = label),position =  position_stack(vjust = 0.5))+
  theme_classic(base_size = 12)+
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5)) + 
  coord_polar(theta = "y", start=0)+
  labs(fill="class",  x=NULL,  y=NULL)+
scale_fill_manual(values=celltype_col)

#fig6h
sigScores <- as.matrix(data.combined@assays$RNA@data)
meta <- data.combined@meta.data[,c("NeoT")]
compare = c("neo-not")
de_gsva <- function(exprSet,meta,compare = NULL,olny_two=T){
  
  allDiff = list()
  design <- model.matrix(~0+factor(meta))
  colnames(design)=levels(factor(meta))
  rownames(design)=colnames(exprSet)
  
  fit <- lmFit(exprSet,design)
  if (olny_two==T) { 
    
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
    
  }else if(length(unique(meta))==21){
    if(is.null(compare)){
      stop("there are 2 Groups,Please set  compare value")
    }
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
    
  }else if(length(unique(meta))>2){
    for(g in colnames(design)){
      fm = ""
      for(gother in colnames(design)[which(!colnames(design) %in% g)]){
        fm = paste0(fm,"+",gother)
      } 
      
      fm = paste0(g,"VsOthers = ",g,"-(",substring(fm,2),")/",ncol(design)-1)
      contrast.matrix <- makeContrasts(contrasts = fm,levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix) 
      fit2 <- eBayes(fit2)
      allDiff[[g]]=topTable(fit2,adjust='fdr',coef=1,number=Inf)
    }
  }else{
    stop("error only have one group")
  }
  
  return(allDiff)
}

Diff =de_gsva(exprSet = sigScores ,meta = meta,compare = compare,olny_two=T)

write.csv(Diff,"Diff_gene_NeoT_vs_not.csv")

VolcanoPlot=function(dif, log2FC=0.2, padj=0.05, 
                     label.symbols=NULL, label.max=90,
                     cols=c("#497aa2", "#ae3137"), title=""){
  if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(dif) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
  }
  rownames(dif)=dif$symbol
  
  # (1) define up and down
  dif$threshold="ns";
  dif[which(dif$log2FoldChange > log2FC & dif$padj <padj),]$threshold="up";
  dif[which(dif$log2FoldChange < (-log2FC) & dif$padj < padj),]$threshold="down";
  dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
  #head(dif)
  #
  tb2=table(dif$threshold); print(tb2)
  library(ggplot2)
  # (2) plot
  g1 = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=0.8, size=2) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    #labs(title= ifelse(""==title, "", paste("DEG:", title)))+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*italic(P.adj)) )+
    ylim(0,80)+
    theme_classic(base_size = 10) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=10)
    ) +
    scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                   paste0("up(",tb2[[3]],')' )),
                       values=c(cols[1], "grey", cols[2]) )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1))); g1;
  # (3)label genes
  if(is.null(label.symbols)){
    dif.sig=dif[which(dif$threshold != "ns" ), ]
    len=nrow(dif.sig)
    if(len<label.max){
      label.symbols=rownames(dif.sig)
    }else{
      dif.sig=dif.sig[order(dif.sig$log2FoldChange), ]
      dif.sig= rbind(dif.sig[1:(label.max/2),], dif.sig[(len-label.max/2):len,])
      label.symbols=rownames(dif.sig)
    }
  }
  dd_text = dif[label.symbols, ]
  print(head(dd_text, n=3))
  print(tail(dd_text, n=3))
  # add text
  library(ggrepel)
  g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                       #size=2.5, 
                       colour="black",alpha=1)
}


diff=Diff[[1]]
dif=data.frame(
  symbol=rownames(diff),
  log2FoldChange=diff$logFC,
  padj=diff$adj.P.Val
)
#dif$padj=dif$padj+1e-153
VolcanoPlot(dif, log2FC=0.2, padj=0.05, title="", label.max = 18,label.symbols=NULL,)

#fig6i
#refer to fig2g









