library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)

#fig10a
source(file.path(code.path, "ML.R"))
FinalModel <- c("panML", "multiLogistic")[1]
Train_expr <- read.csv(file.path(data.path, "Training_expr.csv"), header = T,  row.names = 1,check.names = F,stringsAsFactors = F)
Train_class <- read.csv(file.path(data.path, "Training_class.csv"), header = T, row.names = 1,check.names = F,stringsAsFactors = F)
colnames(Train_class)="outcome"
sam1=rownames(Train_class)[Train_class$outcome==1]
sam2=rownames(Train_class)[Train_class$outcome==0]
Train_class=data.frame(Train_class[c(sam1,sam2),])
rownames(Train_class)=c(sam1,sam2)
colnames(Train_class)="outcome"
comsam <- intersect(rownames(Train_class), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_class <- Train_class[comsam,,drop = F]
Test_expr <- read.csv(file.path(data.path, "Testing_expr.csv"), header = T, row.names = 1,check.names = F,stringsAsFactors = F)
Test_class <- read.csv(file.path(data.path, "Testing_class.csv"), header = T,  row.names = 1,check.names = F,stringsAsFactors = F)
colnames(Test_class)=c("Cohort","outcome")
sam1=rownames(Test_class)[Test_class$outcome==1]
sam2=rownames(Test_class)[Test_class$outcome==0]
Test_class=data.frame(Test_class[c(sam1,sam2),])
rownames(Test_class)=c(sam1,sam2)
colnames(Test_class)=c("Cohort","outcome")
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_class <- Test_class[comsam,,drop = F]
comgene1=rownames(Train_expr)[rowSums(Train_expr)!=0]
comgene2=rownames(Test_expr)[rowSums(Train_expr)!=0]
comgene <- intersect(comgene1,comgene2)
comgene=intersect(comgene,degs$gene[1:50])
Train_expr <- t(Train_expr[comgene,]) 
Test_expr <- t(Test_expr[comgene,]) 
methods <- read.xlsx(file.path(code.path, "methods.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)
classVar = "outcome" 
min.selected.var = 1 
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+") 
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1]) 
preTrain.method = unique(unlist(preTrain.method)) 
preTrain.var <- list()
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, 
                                 Train_set = Train_set, 
                                 Train_label = Train_class, 
                                 mode = "Variable",      
                                 classVar = classVar) 
}
preTrain.var[["simple"]] <- colnames(Train_set)
preTrain.var[["Lasso"]]=degs$gene[1:5]
model <- list() 
Train_set_bk = Train_set
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method 
  method <- strsplit(method, "\\+")[[1]]   
  if (length(method) == 1) method <- c("simple", method)   
  Variable = preTrain.var[[method[1]]] 
  Train_set = Train_set_bk[, Variable]  
  Train_label = Train_class            
  model[[method_name]] <- RunML(method = method[2],      
                                Train_set = Train_set,    
                                Train_label = Train_label,
                                mode = "Model",          
                                classVar = classVar)     
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set = Train_set_bk; rm(Train_set_bk) 
saveRDS(model, file.path(res.path, "model.rds")) 

if (FinalModel == "multiLogistic"){
  logisticmodel <- lapply(model, function(fit){ 
    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family = "binomial", 
               data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit) 
    return(tmp)
  })
}
saveRDS(logisticmodel, file.path(res.path, "logisticmodel.rds")) 

model <- readRDS(file.path(res.path, "model.rds"))
methodsValid <- names(model)
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set)) 
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) 
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set)) 
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
#Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),]) 
write.table(Class_mat, file.path(res.path, "Class_mat.txt"), 
            sep = "\t", row.names = T, col.names = NA, quote = F)
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file.path(res.path, "fea_df.txt"),
            sep = "\t", row.names = F, col.names = T, quote = F)

AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],    
                                Test_set = Test_set,      
                                Test_label = Test_class,  
                                Train_set = Train_set,   
                                Train_label = Train_class, 
                                Train_name = "noot_train",      
                                cohortVar = "Cohort",      
                                classVar = classVar)     
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat, file.path(res.path, "AUC_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)          
avg_AUC <- sort(avg_AUC, decreasing = T)     
AUC_mat <- AUC_mat[names(avg_AUC), ]     
fea_sel <- fea_list[[rownames(AUC_mat)[1]]] 
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3)) 
if(ncol(AUC_mat) < 3) { 
  CohortCol <- c("red","blue") 
} else { 
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") 
}
names(CohortCol) <- colnames(AUC_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, 
                    avg_AUC, 
                    CohortCol, "steelblue", 
                    cellwidth = cellwidth, cellheight = cellheight, 
                    cluster_columns = F, cluster_rows = F)
pdf("fig10a.pdf"), width = 8, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())

#fig10b
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(survival)
library(stringr)
library(viridis)
library(scales)
Coxoutput=data.frame()
for(i in tum_type){
    sur_hall_sub=sur[sur$type==i,]
    realdata=sur_hall_sub
    cox <- coxph(Surv(DSS.time, DSS) ~ sig, data = realdata)
    coxSummary = summary(cox)
    Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                    z=coxSummary$coefficients[,"z"],
                                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"]
}
write.csv(Coxoutput,'cox_output_dss.csv', row.names = F)
#The OS and PFI analysis is the same
library(reshape2)
mycol=tcga_color_system[tumors,]$color
botdf <- cox_output
botdf$group <- factor(botdf$group, levels = tumors)
botdf$dirct <- factor(ifelse(botdf$HR > 1,"Risky","Protective"), levels = c("Risky","Protective"))
botdf$event <- factor(botdf$event, levels = c("OS","PFI","DSS"))
botdf$pvalue=as.numeric(botdf$pvalue)
ggplot(botdf,aes(x = group, y = event,size = -log10(pvalue))) +
    geom_point(shape = 21, aes(size=-log10(pvalue), col = dirct, fill = group), position = position_dodge(0), stroke = 1) +
    scale_fill_manual(values = mycol) +
    scale_color_manual(values = c("grey80","black")) +
    xlab(NULL) + ylab(NULL) +
    labs(size = "-log10(P-value)", col = "HR") +
    scale_size_continuous(range = c(3,8)) + 
    guides(fill = FALSE) + 
    theme_bw() + 
    theme(axis.ticks = element_line(size = 0.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90,size = 13,vjust = 0.5,hjust = 0.5,colour = mycol),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom",
          legend.margin=margin(t= 0, unit='cm'),
          legend.spacing = unit(0.1,"in"))
ggsave("fig10b.pdf",width = 8,height = 4)

#fig10c-e
Neo_sig=Neo_sig[rownames(sur),]
sur$sig=Neo_sig
dat=sur
risk=dat$sig
res.cut=surv_cutpoint(dat,time="OS.time",
                      event ="OS",variables="age")
res.cut=res.cut$cutpoint$cutpoint
risk<-as.vector(ifelse(risk >res.cut,"high","low"))
dat$group<-risk
dat$OS.time=dat$OS.time/365
dat$DSS.time=dat$DSS.time/365
dat$PFI.time=dat$PFI.time/365
dat$DFI.time=dat$DFI.time/365
sur=dat
fitd <- survdiff(Surv(OS.time,OS) ~ group,
                 data      = dat,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS.time, OS)~ group,
               data      = dat,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
ps <- pairwise_survdiff(Surv(OS.time, OS)~ group,
                        data            = dat,
                        p.adjust.method = "none") 

mycol<- c( "#E31A1C","#1F78B4")
names(fit$strata) <- gsub("group=", "", names(fit$strata))
p <- ggsurvplot(fit= fit,
                conf.int FALSE, 
                risk.table= TRUE, 
                risk.table.col= "strata",
                palette= mycol, 
                data= dat,
                size= 1,
                break.time.by= 5, 
                legend.title= "",
                xlab= "Time (years)",
                ylab= "Overall survival",
                risk.table.y.text = FALSE,
                tables.height= 0.3) 
p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 0, y = 0.55,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
#The DSS and PFI analysis is the same

#fig10f
sur$path=sur$sig
x=aggregate(x=sur$path,by=list(type=sur$type),mean)
nets.mean=x$x
names(nets.mean) <- x$type
nets.mean <- sort(nets.mean, decreasing = T) 
tumor.level <- names(nets.mean) 
topdf=sur
topdf$path=topdf$path
mycol=tcga_color_system[tumor.level,]$color
topdf$type <- factor(topdf$type, levels = tumor.level)
topdf$path <- scale(topdf$path) 
ggplot(data = topdf,aes(x = type, y = path, fill = type))+ 
    geom_hline(yintercept = 0, color="black",
               linetype="longdash", lwd = 0.6) +
    geom_boxplot(aes(color=type),fill=NA,width=.8,cex=0.7,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+
    geom_jitter(shape=16, position=position_jitter(0.2),aes(color=type), size = 0.2,alpha=0.8) + 
    geom_violin(aes(color=type), draw_quantiles = c(0.25, 0.5, 0.75),fill=NA, size = 0.6) + 
    stat_summary(fun = "mean", 
                 geom = "point",
                 color = "black") +
    scale_fill_manual(values = alpha(mycol, 0.8)) +   scale_color_manual(values = alpha(mycol, 0.8)) + 
    ylab("Innate") + xlab("") +
    theme_bw() +
    theme(axis.ticks = element_line(size = 0.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 12, color = mycol, angle = 90),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none")
ggsave("fig10f.pdf",width = 8,height = 4)
  
 #fig10g
library(ComplexHeatmap) 
library(circlize) 
library(ChAMPdata) 
library(data.table) 
library(impute)
library(randomcoloR)
library(ggpubr)
library(GSVA)
library(clusterProfiler)
library(dplyr)
data("probe.features")


tumors=c("ACC","THCA","THYM","LGG","GBM",
             "PCPG","UVM","SKCM","SARC","MESO",
             "UCS","UCEC","OV","CESC","BRCA",
             "TGCT","BLCA","PRAD","LIHC","LAML",
             "DLBC","KIRP","KICH","KIRC","ESCA",
             "STAD","PAAD","CHOL","COAD","READ",
             "HNSC","LUSC","LUAD")
frg <- rownames(immunomodulator)
promoter <- probe.features[which(probe.features$feature %in% c("TSS1500","TSS200")),]
promoter <- promoter[which(promoter$gene %in% frg),] 
write.table(promoter, "promoter_annotation_for_interested_genes.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
comgene <- intersect(rownames(expr),frg)
write.csv(comgene,"comm_gene.csv")
rm(expr)
gc()
for (i in tumors) { 
  message ("--",i,"...")
  load(file.path("./meth/", paste0("TCGA-",i,"_methy450.RData"))) 
  methy450 <- as.data.frame(methy450) 
  rownames(methy450) <- methy450[,1] 
  methy450 <- methy450[,-1] 
  dimname <- dimnames(methy450) 
  methy450 <- sapply(methy450, as.numeric) 
  dimnames(methy450) <- dimname 
  methy450 <- as.data.frame(methy450) 
  compb <- intersect(rownames(methy450),rownames(promoter)) 
  methy450 <- methy450[compb,] 
  methy450$gene <- promoter[compb,"gene"] 
  methy450 <- apply(methy450[,setdiff(colnames(methy450), "gene")], 2, function(x) tapply(x, INDEX=factor(methy450$gene), FUN=median, na.rm=TRUE))
  methy450 <- as.data.frame(methy450)
  write.table(methy450, paste0("TCGA_",i,"_methy450_subset.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
  rm(methy450); gc() }
for (clus in c("C1","C2")) {
  meth_for_clus=data.frame(comgene)
  rownames(meth_for_clus)=comgene
  colnames(meth_for_clus)="comgene"
  meth_for_clus=as.data.frame(t(meth_for_clus))
  sam_clus=rownames(sur[sur$group==clus,])
  for (i in tumors) {
    message("--",i,"...")
    meth_subset <- read.table(paste0("TCGA_",i,"_methy450_subset.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    colnames(meth_subset) <- substr(colnames(meth_subset),1,15)
    meth_sam=colnames(meth_subset)
    tumsam_meth <- meth_sam[substr(meth_sam,14,14) == "0"] 
    meth_subset <- meth_subset[,tumsam_meth] 
    colnames(meth_subset) <- substr(colnames(meth_subset),1,12)
    comsam <- intersect(sam_clus,colnames(meth_subset)) 
    if (length(comsam)!= 0) {
      if (length(comsam)== 1) {
         meth_subset2 <- meth_subset[,comsam]
        meth_subset2=as.data.frame(meth_subset2)
        rownames(meth_subset2)=rownames(meth_subset)
        colnames(meth_subset2)=comsam
        meth_subset2$comgene=rownames(meth_subset2)
        meth_subset2=as.data.frame(t(meth_subset2))
        meth_for_clus=bind_rows(meth_for_clus,meth_subset2,.id = "comgene") 
      }else{
        meth_subset <- meth_subset[,comsam]
      meth_subset=as.data.frame(meth_subset)
      colnames(meth_subset)=comsam
      meth_subset$comgene=rownames(meth_subset)
      meth_subset=as.data.frame(t(meth_subset))
      meth_for_clus=bind_rows(meth_for_clus,meth_subset,.id = "comgene")
    }
    }
  }
  write.csv(t(meth_for_clus),paste0("meth_for_",clus,".csv"))
}
meth_for_C1 <- read.csv("./meth_for_C1.csv", row.names=1)
meth_for_C2 <- read.csv("./meth_for_C2.csv", row.names=1)
meth=cbind(meth_for_C1,meth_for_C2)
corExpMeth <- 
  as.data.frame(matrix(NA,
                       nrow = length(comgene),
                       ncol = 5, 
                       dimnames = list(comgene, 
                                       unique(c("C1","C2")))))
expr_sub2=expr_sub
colnames(expr_sub2)= gsub("-",".",colnames(expr_sub2))

for (i in comgene) {
  if(!is.element(i, rownames(expr_sub2)) | !is.element(i, rownames(meth))) {
    corExpMeth[i,] <- NA 
  } else { 
    for (j in c("C1","C2")) {
      sam_all <- rownames(sur[which(sur$group == j),,drop = F])
      sam_all=gsub("-",".",sam_all)
      comsam=intersect(sam_all,colnames(meth))
      expr.subset <- as.numeric(expr_sub2[i, comsam])
      meth.subset <- as.numeric(meth[i, comsam])
      if (sum(is.na(meth.subset))!=length(meth.subset)) {
        ct <- cor.test(expr.subset, meth.subset, method = "spearman",na.action = "na.exclude" ) 
        corExpMeth[i, j] <- ct$estimate
      }else{ corExpMeth[i, j] = NA}
    }
  }
}
write.csv(corExpMeth,"corExpMeth.csv")
expMat <- as.data.frame(t(expr_sub))
sur=sur[colnames(expr_sub),]
expMat$group <- sur$group
expMat <- as.data.frame(t(apply(expMat[,setdiff(colnames(expMat), "group")], 2, 
                                function(x) 
                                  tapply(x, 
                                         INDEX = factor(expMat$group), 
                                         FUN = mean, 
                                         na.rm = TRUE)))) 
write.table(expMat,"expMat.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
cna <- fread("./Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")
cna <- as.data.frame(cna); rownames(cna) <- cna[,1]; cna <- cna[,-1]
cna$gene <- sapply(strsplit(rownames(cna),"|",fixed = T),"[",1)
cna <- cna[!duplicated(cna$gene),]; cna <- cna[,setdiff(colnames(cna),"gene")]
is.element(comgene,rownames(cna)) 
cna <- cna[intersect(rownames(cna),comgene),]
cna[cna > 1] <- 1 
cna[cna < -1] <- -1
colnames(cna)=substr(colnames(cna),1,12)
comsam <- intersect(colnames(expr_sub), colnames(cna))
cna <- cna[,comsam]
write.table(cna, file = "cna.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
ampFreq <- delFreq <- 
  as.data.frame(matrix(NA,
                       nrow = length(comgene),
                       ncol = 5, 
                       dimnames = list(comgene,unique(c("C1","C2")))))
for (i in comgene) {
  if(!is.element(i, rownames(cna))) { 
    ampFreq[i,] <- NA 
    delFreq[i,] <- NA
  } else { 
    ampFreqInAll <- sum(as.numeric(cna[i,]) == 1)/ncol(cna) 
    delFreqInAll <- sum(as.numeric(cna[i,]) == -1)/ncol(cna) 
    for (j in c("C1","C2")) {
      sam <- rownames(sur[which(sur$group == j),,drop = F])
      comsam=intersect(sam,colnames(cna))
      cna.subset <- cna[, comsam]
      ampFreqInSubt <- sum(as.numeric(cna.subset[i,]) == 1)/length(comsam) 
      delFreqInSubt <- sum(as.numeric(cna.subset[i,]) == -1)/length(comsam) 
      ampFreqInDiff <- ampFreqInSubt - ampFreqInAll 
      delFreqInDiff <- delFreqInSubt - delFreqInAll 
      ampFreq[i, j] <- ampFreqInDiff
      delFreq[i, j] <- delFreqInDiff
    }
  }
}
write.table(ampFreq,"ampFreq.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(delFreq,"delFreq.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
library(ComplexHeatmap) 
library(circlize)
annCol <- data.frame(subtype = c("C1","C2"),
                     row.names = c("C1","C2"))
annColors <- list()
annColors[["subtype"]] <- c("C1" = "#1F78B4",
                            "C2" = "#E31A1C",
                            "C1" = "#33A02C",
                            "C5" = "#6A3D9A",
                            "C3" = "#E31A1C")
top_anno <- HeatmapAnnotation(df= annCol,
                              col= annColors,
                              gp= gpar(col = "grey80"), 
                              simple_anno_size= unit(3.5, "mm"), 
                              show_legend= F, 
                              show_annotation_name = F, 
                              border= FALSE) 
annRow <- immunomodulator[comgene,]
annRow[which(annRow$Category == "Co-stimulator"),"Category"] <- "Co-stm" 
annRow[which(annRow$Category == "Co-inhibitor"),"Category"] <- "Co-ihb"
annRow[which(annRow$Category == "Cell adhesion"),"Category"] <- "Cell\nadhesion" 
annRow[which(annRow$Category == "Antigen presentation"),"Category"] <- "Antigen\npresentation"
annRow$Category <- factor(annRow$Category, levels = c("Co-stm","Co-ihb","Ligand","Receptor","Cell\nadhesion","Antigen\npresentation","Other")) 
annRow$ICI <- factor(annRow$ICI, levels = c("Inhibitory","N/A","Stimulatory"))
annRowColors <- list("ICI" = c("Inhibitory" = "black","N/A" = "#888888","Stimulatory" = "#E59E02"))
left_anno <- HeatmapAnnotation(df= annRow[,"ICI",drop = F],
                               which= "row", 
                               gp= gpar(col = "grey80"), 
                               col= annRowColors,
                               simple_anno_size= unit(3.5, "mm"), 
                               show_annotation_name = F,
                               border= F)
left_anno <- HeatmapAnnotation(df= annRow[,"cris",drop = F],
                               which= "row",
                               gp= gpar(col = "grey80"), 
                               simple_anno_size     = unit(3.5, "mm"), 
                               show_annotation_name = F,
                               border= F)
draw(left_anno)
heatmap.BlWtRd <- c("#1F66AC", "#75AFD3", "grey90", "#FAB99B", "#B2192B")
exp <- apply(expMat, 1, scale)
rownames(exp) <- colnames(expMat)
exp <- t(exp)
exp=exp[,rownames(annCol)]
col_expr <- colorRamp2(seq(min(na.omit(exp)), max(na.omit(exp)), length = 5), heatmap.BlWtRd)
hm.expr <- Heatmap(matrix= as.matrix(exp),
                   col= col_expr,
                   border= NA, 
                   rect_gp = gpar(col = "grey80"), 
                   cluster_rows= F, 
                   cluster_columns= F, 
                   show_row_names= T, 
                   row_names_side= "left", 
                   row_names_gp= gpar(fontsize = 10), 
                   show_column_names  = F, 
                   column_names_side  = "top", 
                   row_split= annRow$Category, 
                   top_annotation= top_anno,
                   left_annotation= left_anno, 
                   name= "mRNA\nExpression", 
                   width= ncol(expMat) * unit(4, "mm"), 
                   height= nrow(expMat) * unit(3.5, "mm")) 
corExpMeth=corExpMeth[,rownames(annCol)]
col_corExprMeth <- colorRamp2(seq(min(na.omit(corExpMeth)), max(na.omit(corExpMeth)), length = 5), heatmap.BlWtRd)
hm.corExprMeth <- Heatmap(matrix= as.matrix(corExpMeth),
                          col= col_corExprMeth,
                          border= NA,
                          rect_gp = gpar(col = "grey80"),
                          cluster_rows= F,
                          cluster_columns= F,
                          show_row_names= F,
                          row_names_side= "left",
                          row_names_gp= gpar(fontsize = 10),
                          show_column_names  = F,
                          column_names_side  = "top",
                          row_split= annRow$Category,
                          row_title= NULL,
                          top_annotation= top_anno,
                          name= "Expression\nvs. Methylation",
                          width= ncol(expMat) * unit(4, "mm"),
                          height= nrow(expMat) * unit(3.5, "mm"))
ampFreq=ampFreq[,rownames(annCol)]
col_ampFreq <- colorRamp2(seq(min(na.omit(ampFreq)), max(na.omit(ampFreq)), length = 5), heatmap.BlWtRd)
hm.ampFreq <- Heatmap(matrix= as.matrix(ampFreq),
                      col= col_ampFreq,
                      border= NA,
                      rect_gp = gpar(col = "grey80"),
                      cluster_rows= F,
                      cluster_columns= F,
                      show_row_names= F,
                      row_names_side= "left",
                      row_names_gp= gpar(fontsize = 10),
                      show_column_names  = F,
                      column_names_side  = "top",
                      row_split= annRow$Category,
                      row_title= NULL,
                      top_annotation= top_anno,
                      name= "Amplification\nFrequency",
                      width= ncol(expMat) * unit(4, "mm"),
                      height= nrow(expMat) * unit(3.5, "mm"))
delFreq=delFreq[,rownames(annCol)]
col_delFreq <- colorRamp2(seq(min(na.omit(delFreq)), max(na.omit(delFreq)), length = 5), heatmap.BlWtRd)
hm.delFreq <- Heatmap(matrix= as.matrix(delFreq),
                      col= col_delFreq,
                      border= NA,
                      rect_gp = gpar(col = "grey70"),
                      cluster_rows= F,
                      cluster_columns= F,
                      show_row_names= F,
                      row_names_side= "left",
                      row_names_gp= gpar(fontsize = 10),
                      show_column_names  = F,
                      column_names_side  = "top",
                      row_split= annRow$Category,
                      row_title= NULL,
                      top_annotation= top_anno,
                      name= "Deletion\nFrequency",
                      width= ncol(expMat) * unit(4, "mm"),
                      height= nrow(expMat) * unit(3.5, "mm"))
pdf(file = "fig10g.pdf", width = 8,height = 12)
draw(hm.expr + hm.corExprMeth + hm.ampFreq + hm.delFreq, 
     heatmap_legend_side = "bottom")
invisible(dev.off())

#fig10h
library(stats)
library(Hmisc)
library(vegan)  library(ggcor) 
library(ggplot2)
library(dplyr)
xcell=ss_imm
data.new <-as.data.frame(t(xcell))  
data.new$sig=sur[rownames(data.new),]$sig
cor_fin=data.frame()
for (i in c(1:28)) {
    cor=cor.test(data.new[,i],data.new[,29])
    cor_fin[i,1]=cor$estimate
    cor_fin[i,2]=cor$p.value
}
rownames(cor_fin)=colnames(data.new)[1:28]
colnames(cor_fin)=c("r","p.value")
cor_fin$spec="sig"
cor_fin$env=rownames(cor_fin)

r.p.data.plot <- cor_fin %>% 
    mutate(r.sign = cut(r, breaks = c(-Inf, 0, Inf), 
                        labels = c("Negative", "Positive")),
           p.sign = cut(p.value, breaks = c(0, 0.05, Inf), 
                        labels = c("P<0.05", "P>=0.05"),
                        include.lowest = TRUE,
                        right = FALSE), 
           r.abs = cut(abs(r), breaks = c(0, 0.1, 0.3, 0.5,1),
                       labels = c("<0.1","0.1-0.3", "0.3-0.5","0.5-1"),
                       include.lowest = TRUE,
                       right = FALSE) )  
r=r.p.data.plot$r
p.value=r.p.data.plot$p.value
r.p.data.plot=r.p.data.plot[,-1]
r.p.data.plot=r.p.data.plot[,-1]
r.p.data.plot$r=r
r.p.data.plot$p.value=p.value
quickcor(data.new, type = "upper", show.diag = F ) + 
    geom_square() +  
    anno_link(data = r.p.data.plot,   
              aes(colour = r.sign,
                  size = r.abs,
                  linetype = p.sign),width = 3, 
              nudge_x = 0,
              curvature = 0.2,alpha=0.4) + 
    scale_size_manual(values = c("<0.1" = 0.25,
                                 "0.1-0.3" = 0.5,
                                 "0.3-0.5" = 1,
                                 "0.5-1" = 2)) +  
    scale_colour_manual(values = c("Negative" = "#1F78B4",
                                   "Positive" = "#B2192B"))+  
    scale_linetype_manual(values = c("P<0.05" = "solid",
                                     "P>=0.05" = "dashed"))+
    scale_fill_gradientn(colours = rev(c("#5b0018", "#e0745a", "#fbfbf7", "#63a1cb", "#052452")),  
                         breaks = seq(-1, 1, 0.2),
                         limits = c(-1, 1))+ 
    guides(
        fill = guide_colorbar(title = "Pearson's r"), 
        linetype = guide_legend(title = NULL),
        colour = guide_legend(title = NULL),
        size = guide_legend(title = "|LMMs r|") 
    )+
    theme(
        legend.key.size = unit(3.5, "mm"), 
        legend.spacing = unit(5, "mm"),
        legend.key = element_blank()     )

#fig10i
ss_imm=as.data.frame(t(ss_imm))
ss_imm$sig=sur$sig
corCnaExpr <- NULL
for (i in tumors) {
  message("--",i,"...")
   sam=rownames(sur[sur$type==i,])
  ss_imm_sub=ss_imm[sam,]
  ss_imm_sub2=ss_imm_sub[,-length(colnames(ss_imm_sub))]
  corTab <- NULL
  for (j in colnames(ss_imm_sub2)) {
    tmp1 <- as.numeric(ss_imm_sub2[,j]) 
    tmp2 <- as.numeric(ss_imm_sub[,length(colnames(ss_imm_sub))]) 
    cor.res <- cor.test(tmp1,tmp2, method = "spearman") 
    corTab <- rbind.data.frame(corTab,
                               data.frame(gene = j,
                                          tumor = i,
                                          Correlation = ifelse(is.na(cor.res$estimate), 0, cor.res$estimate),
                                          Pvalue = ifelse(is.na(cor.res$p.value), 1, cor.res$p.value),
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
  }
  corCnaExpr <- rbind.data.frame(corCnaExpr,
                                 corTab,
                                 stringsAsFactors = F)
}
ggplot(corCnaExpr, aes(x=tumor,y=gene)) +
  geom_point(aes(size=-log10(Pvalue),color=Correlation)) +
  scale_color_gradientn('Correlation',colours = c('#4575B4','#8DBBD9','#FBFDC8','#EF6E46','#D73027'))+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, size = 12, hjust = 0.9, vjust = 0.9, color = col_tum),
     axis.title = element_blank(),
    panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
    legend.position = "right",
    plot.margin = unit(c(1,1,1,1), "lines"))















































