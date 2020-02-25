
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq")
#biocLite("limma")
#install.packages("gplots")



setwd("C:\\Users\\think\\Desktop\\RNA")           
library("DESeq")
library("limma")

rt=read.table("WL1A8_count.txt",sep="\t",header=T)  
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=round(data,0)

#View(head(data))

group=c(rep("WT",3),rep("KO",3))  
design = factor(group)
newTab = newCountDataSet( data, design )
newTab = estimateSizeFactors(newTab)
newData=counts(newTab, normalized=TRUE )

#have replicates
newTab = estimateDispersions( newTab, fitType = "local")
diff = nbinomTest( newTab, "WT", "KO")
diff = diff[is.na(diff$padj)==FALSE,]
diff = diff[order(diff$pval),]
write.table( diff, file="DESeqOut.xls",sep="\t",quote=F,row.names=F)


normalizeExp=rbind(id=colnames(newData),newData)
diffExp=rbind(id=colnames(newData),newData[diff$id,])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         
