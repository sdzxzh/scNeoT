library(Seurat)
library(scRepertoire)

#fig5a
a_NS <- read.csv("./data/NS/filtered_contig_annotations.csv")
b_DC <- read.csv("./data/DC/filtered_contig_annotations.csv")
c_OCDC <- read.csv("./data/OCDC/filtered_contig_annotations.csv")
d_T <- read.csv("./data/T/filtered_contig_annotations.csv")
e_NeOT <- read.csv("./data/NeOT/filtered_contig_annotations.csv")
f_Noot <- read.csv("./data/Noot/filtered_contig_annotations.csv")
contig_list <- list(a_NS,b_DC,c_OCDC,d_T,e_NeOT)
combined <- combineTCR(contig_list, 
                        samples = c("a_NS", "b_DC", "c_OCDC", "d_T","e_NeOT"), 
                        ID =  c("a_NS", "b_DC", "c_OCDC", "d_T","e_NeOT"), 
                        cells ="T-AB")
pdf("fig5a.pdf",8,3)
clonalDiversity(combined, cloneCall = "aa", group = "ID") +scale_fill_manual(values = col)
dev.off()

#fig5b
seurat1 <- combineExpression(combined, data.combined, cloneCall="aa")
pdf("fig5b.pdf",5,5)
alluvialClonotypes(seurat1, cloneCall = "aa", 
                   y.axes = c("group", "fincell"), 
                   color = "celltype") +scale_fill_manual(values =celltype_col )
dev.off()

#fig5c
subset <- subset(seurat1, group == "group_name")
circles <- getCirclize(subset, group.by = "celltype")
grid.cols <- scales::hue_pal()(length(unique(subset@active.ident)))
names(grid.cols) <- levels(subset@active.ident)
pdf("fig5c.pdf",7,7)
chordDiagram(circles, self.link = 1, 
             grid.col = celltype_col, directional = 1, 
             direction.type =  "arrows",
             link.arr.type = "big.arrow")
dev.off()

#fig5f-j
#Prepare network files
#The input data is the output of GLIPH2 analysis
library(foreach)
library(doParallel)
cl <- makeCluster(14)
registerDoParallel(cl)
data=input
index=unique(data$index)
length=length(index)-1
fin<-foreach(i =  1:length,.combine = rbind) %dopar% {
  index_source=index[i]
  index_target=index[i+1:length(index)]
  index_target=na.omit(index_target)
  row=0
  res=data.frame()
  for (target in c(1:length(index_target))) {
    row=row+1
    count=length(intersect(data$TcRb[data$index==index_source],
                           data$TcRb[data$index==index_target[target]]))
    res[row,1]=index_source
    res[row,2]=index_target[target]
    res[row,3]=count
  }
  return(res)
}
stopCluster(cl)
res_sub=fin[fin$V3!=0,]
colnames(res_sub)=c("Source","Target","Weight")
write.csv(res_sub,"res.csv")
#Then visualize the network using Gephi software

#fig5k
#rowann is the row annotation of the index
#colann is the row annotation of the CDR3s
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)
library(Biostrings)
library(motifStack)
library(MotifDb)
dat=dcast(long_dat, index ~ cdr3)
row.names(dat)=dat[,1]
dat=dat[,-1]
onco.input=dat
onco.input[onco.input == 1] <- "Mutated" 
onco.input[onco.input == 2] <- "Mutated2"
onco.input[onco.input == 3] <- "Mutated3" 
onco.input[onco.input == 0] <- "" 
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
  },
  Mutated = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#FA9A99", col = "#FA9A99")) 
  },
  
  Mutated2 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#E21A1C", col = "#E21A1C")) 
  },
  
  Mutated3 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#A60000", col = "#A60000")) 
  }
)
col = c("Mutated" ="#FA9A99",
        "Mutated2" ="#E21A1C",
        "Mutated3" ="#A60000") 
group=rowann$group
ifng_noot = rowann$ifng_noot
ifng_other = rowann$ifng_other
clone = rowann$clone
names(group)=names(ifng_noot)=names(ifng_other)=names(clone)<-rownames(rowann)
left_ann=rowAnnotation(group = group,
                          ifng_noot = ifng_noot,
                          ifng_other = ifng_other,
                          clone = clone,
                          col = list(group = c("Mix" = "grey30",
                                               "Noot" = "grey90",
                                               "NS" = "#A48AD3",
                                               "Dc" = "#1CC5FE",
                                               "OCDC" = "#6FC7CF",
                                               "T" = "#FBA27D",
                                               "NeoT" = "#FB7D80"),
                                     ifng_noot = c("Y" = "#F768A1",
                                                   "N" = "grey90"),
                                     ifng_other= c("Y" = "#F768A1",
                                                   "N" = "grey90"),
                                     clone= c("Clone" = "#EB7C20",
                                                   "non_Clone" = "grey90")
                          ),
                          width = ncol(rowann_sub)*unit(0.1, "mm")
                        )
top_ann <- HeatmapAnnotation(panpep=anno_barplot(as.matrix(colann),
                                                border = FALSE,
                                                gp = gpar(fill = c(brewer.pal(6,"Paired")[c(2,1,6,5)],"grey90"), 
                                                          border = NA, 
                                                          lty = "blank"),
                                                height = unit(3, "cm")) )
draw(top_ann)
order=names(sort(rowSums(colann),decreasing = T)) 
heat=oncoPrint(onco.input[,order], 
          alter_fun = alter_fun,show_pct = F,top_annotation = top_ann,
          right_annotation = right_ann)
pdf("fig5k_1.pdf", height = 7 , width =9)
print(heat)
dev.off()
#neo_trb is the CDR3s of neoantigen-reactive T
dat=neo_trb
colnames(dat)=c("index","TcRb")
motifs=list()
motifs_pcm=list()
id=1
for (i in unique(dat$index)) {
  dat_sub=dat[dat$index==i,]
  aa=strsplit(dat_sub$TcRb,"")
  aa2=data.frame(aa)
  aa3=unique(unlist(aa))
  lenth_aa=length(aa[[1]])
  aa2=as.data.frame(t(aa2))
  tmp_list=list()
  for (t in c(1:lenth_aa)) {
    position_sub=aa2[,t]
    stat=table(position_sub)
    stat=as.data.frame(stat)
    tmp=c()
    for (a in aa3) {
      if (stat$position_sub==a) {
        tmp=c(tmp,stat[stat$position_sub==a,]$Freq)
      }else{tmp=c(tmp,0)}
    }
    tmp_list[[t]]=tmp
  }
  motif=data.frame(tmp_list)
  rownames(motif)=aa3
  colnames(motif)=c(1:lenth_aa)
  motif<-pcm2pfm(motif)
  motif_pfm<-new("pfm", mat=motif, name=paste0("id",i), 
             color=colorset(alphabet="AA",colorScheme="chemistry"))
  motifs[[id]]=motif_pfm
  motif_pcm<-new("pcm", mat=as.matrix(motif), name=paste0("id",i),
                 color=colorset(alphabet="AA",colorScheme="chemistry"))
  motifs_pcm[[id]]=motif_pcm
  id=id+1
}
pfms=motifs
#the logo plot
plot(pfms[[the index we need]])















