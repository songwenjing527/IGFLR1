library("dplyr")
library("DESeq2")
library("limma")
library("edgeR")

rt=read.table("symbol1.txt",sep="\t",header=T,check.names=F) 
#rt %>% distinct(id, .keep_all = TRUE)
index <- duplicated(rt[,1])
index
rt1 <- rt[!index,]
View(rt1)
rownames(rt1)=rt1[,1]
exp=rt1[,2:ncol(rt1)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data[data<0]=0
View(data)

table(substr(colnames(data),14,14))  
t_index <- which(substr(colnames(data),14,14) == '0')
n_index <- which(substr(colnames(data),14,14) == '1')
data <- as.data.frame(data)

exp_tumor <- select(data,t_index)
View(exp_tumor)
write.csv(exp_tumor,file="exp_tumor.csv")
exp_t <- t(exp_tumor)
View(exp_t)
exp_normal <- select(data,n_index)
View(exp_normal)
exp_tn <- t(exp_normal)
View(exp_tn)
exp_t <- as.data.frame(exp_t)
exp_tn <- as.data.frame(exp_tn)
t <- exp_t["IGFLR1"]
n <- exp_tn["IGFLR1"]
View(t)
View(n)
w <- rbind(t, n)
View(w)
write.csv(w,file="exp_IGFLR1.csv")
exp_t1 <- as.data.frame(exp_t)
exp_t2 <- exp_t1[order(exp_t1$IGFLR1),]
exp_t3 <- t(exp_t2)
View(exp_t3)
#median(exp_t1$IGFLR1,na.rm = T)=0.8907632

#group
group_list <- c(rep('low',268),rep('high',267))
group_list <- factor(group_list)
design <- model.matrix(~0+group_list)
design
rownames(design) = colnames(exp_t3)
colnames(design) <- levels(group_list)
exp_t3 <- as.matrix(exp_t3)
#Differential expression matrix
DGElist <- DGEList( counts = exp_t3, group = group_list )
keep_gene <- rowSums( cpm(DGElist) > 0 ) >= 1 # 自定义
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]

DGElist <- calcNormFactors( DGElist )
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(contrasts = c('low-high'), levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

nrDEG_limma_voom = topTable(fit2, coef = 'low-high', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
View(nrDEG_limma_voom)
padj = 0.05 
foldChange= 1 

nrDEG_limma_voom$sig[nrDEG_limma_voom$P.Value > padj|nrDEG_limma_voom$logFC < foldChange| nrDEG_limma_voom$logFC > -foldChange] <- "no"
nrDEG_limma_voom$sig[nrDEG_limma_voom$P.Value <= padj & nrDEG_limma_voom$logFC >= foldChange] <- "down"
nrDEG_limma_voom$sig[nrDEG_limma_voom$P.Value <= padj & nrDEG_limma_voom$logFC <= -foldChange] <- "up"
View(nrDEG_limma_voom)
write.csv(nrDEG_limma_voom, file = 'limma_all.csv')

nrDEG_limma_voom <- read.csv(file="limma_all1.csv",sep = ',',header=T,row.names=1)
nrDEG_limma_voom_signif = nrDEG_limma_voom[(nrDEG_limma_voom$adj.P.Val < padj & 
                                              (nrDEG_limma_voom$logFC1>foldChange | nrDEG_limma_voom$logFC1<(-foldChange))),]
nrDEG_limma_voom_signif = nrDEG_limma_voom_signif[order(nrDEG_limma_voom_signif$logFC1),]
write.csv(nrDEG_limma_voom_signif, file = 'limma_signif.csv')

outDiff <- read.csv(file="limma_signif.csv",sep = ',',header=T)
exp <-  read.csv(file="exp_tumor.csv",sep = ',',header=T)
View(outDiff)
View(exp)
x <- semi_join(exp,outDiff,by="gene")
View(x)
write.csv(x,file="DiffbExp.csv",sep=",",col.names=F,quote=F)

#volcano plot
nrDEG_limma_voom <- read.csv(file="limma_all1.csv",sep = ',',header=T,row.names=1)
library(ggplot2)
p<-ggplot(nrDEG_limma_voom,aes(x=nrDEG_limma_voom$logFC1,y=-log10(nrDEG_limma_voom$P.Value),colour=sig))+xlab("log2 Fold Change")+ylab("-log10P-Value")+
  geom_point(size=5,alpha=0.6,cex.axis=5, cex.lab=5)+
  theme(axis.title.x = element_text(size=25),
        axis.text.x=element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=25),
        legend.title=element_text(family="Times",size=25),
        legend.text=element_text(family="Times",size=25))+
  scale_color_manual(values =c("#0072B5","black","#BC3C28"))+#color
  geom_vline(aes(xintercept=1), colour="gray",size=1.2 ,linetype=2)+
  geom_vline(aes(xintercept=-1), colour="gray",size=1.2 ,linetype=2)+
  geom_hline(aes(yintercept=-log10(0.05)),colour="gray",size=1.2 ,linetype=2) 
guides(fill=guide_legend(title=NULL))
print(p)
ggsave(("limma_vol.tiff"),plot = print(p),width = 10,height = 8,units = "in")

#heatmap
outDiff <- read.csv(file="DiffbExp.csv",sep = ',',header=T,row.names=1)
pdf(file="heatmap.pdf",height=12,width=20)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = F,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=10)
dev.off()

#GO and KEGG
#ID change
library("org.Hs.eg.db")
rt=read.table("id.txt",sep="\t",check.names=F,header=T)
View(rt)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id1.txt",sep="\t",quote=F,row.names=F)

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

#GO
rt=read.table("id1.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID

kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)

pdf(file="barplotGO.pdf",width = 15,height = 7)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

pdf(file="bubbleGO.pdf",width = 14,height = 7)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#KEGG
rt=read.table("id1.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)

pdf(file="barplotKEGG.pdf",width = 9,height = 6)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()

pdf(file="bubbleKEGG.pdf",width = 8,height = 6)
dotplot(kk, showCategory = 30)
dev.off()

#GSEA
#GSEA input file
c <- read.csv(file="Intersection gene.csv",sep = ',',header=T)
View(c)

e <- read.csv(file="exp_tumor.csv",sep = ',',header=T)
View(e)
as.character(c$gene)
as.character(e$gene)
i <- semi_join(e,c,by="gene")
View(i)
rownames(i)=i[,1]
i=i[,-1]
t <- t(i)
t <- as.data.frame(t)
View(t)
median(t$IGFLR1,na.rm = T)
a <- ifelse(t$IGFLR1>median(t$IGFLR1,na.rm = T),"high","low")
t$IGFLR1 <- a
View(t)
write.csv(i,file="Intersection gene expression.csv",row.names=F)
x <- t(t)
write.csv(x,file="Intersection gene expression1.csv",row.names=T)

#multiple pathway

library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

files=grep(".xls",dir(),value=T)                                        
data = lapply(files,read.delim)                                         
names(data) = files

dataSet = ldply(data, data.frame)
dataSet$pathway = gsub(".xls","",dataSet$.id)                           

gseaCol=c("#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,fill=pathway,group=pathway))+
  geom_point(shape=21) + scale_fill_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits =c(min(dataSet$RUNNING.ES-0.02), max(dataSet$RUNNING.ES+0.02))) +   
  theme_bw() + theme(panel.grid =element_blank()) + theme(panel.border = element_blank()) + 
  theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(), axis.title.y = element_text(size=15),axis.text.y = element_text(size=15),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) + guides(fill=guide_legend(title = NULL)) + 
  theme(legend.background = element_blank(),legend.title=element_text(size=10),legend.text=element_text(size=10)) + theme(legend.key = element_blank())
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "High Exp<----------->Low Exp", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.x = element_text(size=15),axis.text.x = element_text(size=15),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)

gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()

pdf('im.pdf',      
    width=11,                
    height=5)               
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
dev.off()

library(limma)
cli <- read.csv("cli_tumor1.csv",header=T,check.names=F)
gene <- read.csv("IGFLR1.csv",header=T,check.names=F)

m <- merge(cli,gene,by="id")
write.csv(m,"cox.csv",row.names = F)

#uniCOX
install.packages('survival')
library(survival)
pFilter=0.05                                                      
setwd("D:/临床+生信/IGFLR1免疫分析/单基因分析修-swj/单因素多因素分析")         
rt=read.csv("cox.csv",header=T,check.names=F,row.names=1)     

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outTab,file="UniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="UniSigExp.txt",sep="\t",row.names=F,quote=F)

rt <- read.table("UniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

pdf(file="uniforest.pdf", width = 7,height =4.5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(5,3))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=1,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#multivariate Cox analysis

library(survival)                                        
rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)    
rt$futime=rt$futime/365

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="Risk1.txt",
            sep="\t",
            quote=F,
            row.names=F)

rt <- read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")

Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

pdf(file="multiforest.pdf", width = 7,height =4.5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(5,3))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=1,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

##nomogram
library(rms)
rt=read.table("Risk1.txt",header=T,sep="\t",row.names=1) 
rt <- rt[,-(9:10)]
dd <- datadist(rt)
options(datadist="dd")
f <- cph(Surv(futime, fustat) ~ age+cancer_status+pathological_M+platelet+IGFLR1, x=T, y=T, surv=T, data=rt, time.inc=3)
surv <- Survival(f)

nom <- nomogram(f, fun=list(function(x)1/(1+exp(-x)),function(x) surv(3, x), function(x) surv(5, x), function(x) surv(10, x)), 
                lp=F, funlabel=c("Risk","3-year survival", "5-year survival", "10-year survival"), 
                maxscale=100, 
                fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
pdf(file="Nomogram.pdf",height=7.5,width=11)
plot(nom)
dev.off()

#C=0.194
validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(futime, fustat) ~ predict(f), data = rt)

#calibration curve
#3-year survival calibration curve
f3<- cph(Surv(futime, fustat) ~ age+cancer_status+pathological_M+platelet+IGFLR1, x=T, y=T, surv=T, data=rt, time.inc=3)
cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3, m=53, B=1000)
pdf(file="cal3.pdf",height=5,width=5)
plot(cal3)
dev.off()

#5-year survival calibration curve
f5<- cph(Surv(futime, fustat) ~ age+cancer_status+pathological_M+platelet+IGFLR1, x=T, y=T, surv=T, data=rt, time.inc=5)
cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5, m=53, B=1000)
pdf(file="cal5.pdf",height=5,width=5)
plot(cal5)
dev.off()

#10-year survival calibration curve
f10<- cph(Surv(futime, fustat) ~ age+cancer_status+pathological_M+platelet+IGFLR1, x=T, y=T, surv=T, data=rt, time.inc=10)
cal10 <- calibrate(f10, cmethod="KM", method="boot", u=10, m=53, B=1000)
pdf(file="cal10.pdf",height=5,width=5)
plot(cal10)
dev.off()

#lasso
install.packages("glmnet")
library("glmnet")
library("survival")

setwd("D:/临床+生信/IGFLR1免疫分析/单基因分析修-swj/单因素多因素分析")    
rt=read.table("UniSigExp.txt",header=T,sep="\t",row.names=1)          
rt$futime=rt$futime/365

rt <- subset(rt, futime >0)
x=as.matrix(rt[,c(3:ncol(rt))])

y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000)

coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="Risk.txt",sep="\t",quote=F,row.names=F)

#survival
library(survival)
library("survminer")

bioSurvival=function(inputFile=null,outFile=null){
  rt=read.table(inputFile,header=T,sep="\t")               
  
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime,fustat) ~ risk, data = rt)
  
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=TRUE,
                     pval=paste0("p=",pValue),
                     pval.size=4,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("red", "blue"),
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="Risk1.txt",outFile="Risk1.pdf")

#risk plot
install.packages("pheatmap")
library(pheatmap)
library(survival)
library("survminer")

rt=read.table("Risk1.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

pdf(file="survival1.pdf",onefile = FALSE,
    width = 5.5,            
    height =5)             
ggsurvplot(fit, 
           data=rt,
           conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()

summary(fit)  

#ROC
install.packages("survivalROC")
library(survivalROC)

rt=read.table("Risk1.txt",header=T,sep="\t",check.names=F,row.names=1)    
pdf(file="ROC1.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",sprintf("%.3f",roc$AUC),")"),
     lwd = 3, cex.main=1.5, cex.lab=1.5, cex.axis=1.5, font=1.5)
abline(0,1)
dev.off()

#risk plot
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file="riskScore1.pdf",width = 10,height = 4)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
dev.off()

#survial status plot
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStat1.pdf",width = 10,height = 4)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

#risk heatmap
install.packages("pheatmap")
library(pheatmap)
rt=read.table("Risk1.txt",sep="\t",header=T,row.names=1,check.names=F)       
rt=rt[order(rt$riskScore),]  
rt1=rt[c(3:(ncol(rt)-2))]
rt1=log2(t(rt1)+0.1)
annotation=data.frame(risk=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmaprisk1.pdf",width = 10,height = 4)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         color = colorRampPalette(c("green", "black", "red"))(50) )
dev.off()
