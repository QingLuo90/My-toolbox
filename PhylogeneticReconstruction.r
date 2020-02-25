library(ape)
library(ggtree)
setwd("D:\\0530\\SciClone\\M1")
k = read.table("k.txt",header=T,row.names = 1)
d = nj(dist.gene(t(k),pairwise.deletion = F))
plot.phylo(d,edge.color="darkblue",edge.width = 1)



###############################################################
library(ape)
library(ggtree)

setwd("E:\\Simon\\Mycode\\KI_0603\\GSE117542_RAW")
rt = read.table("PGC_AlleleTable_for_tree.txt",header=T)
#table(rt$Cell)
#E1 E2 E3 E5 E6 E7 
#70 88 54 56 96 49 

E1 = rt[which(rt[,1]=="E1"),]
E1 = E1[,-1]
rownames(E1)=E1[,1]
E1 = E1[,-1]

sort(table(E1$allele1),decreasing = T)

E1_1 = E1[which(E1$allele1=="139:11D"),]
E1_2 = E1[which(E1$allele1=="104:51D"),]
E1_3 = E1[which(E1$allele1=="47:92D"),]
E1_4 = E1[which(E1$allele1=="136:1D"),]


k1=t(E1_1)
d1 = nj(dist.gene(t(k1),pairwise.deletion = F))
ggtree(d1 ,ladderize = T)+ geom_point(col="royalblue")+ theme(legend.position = "right")+geom_tiplab(size=2,col="royalblue")

k1=t(E1_2)
d1 = nj(dist.gene(t(k1),pairwise.deletion = F))
ggtree(d1 ,ladderize = T)+ geom_point(col="#771111")+ theme(legend.position = "right")+geom_tiplab(size=2,col="#771111")

k1=t(E1_3)
d1 = nj(dist.gene(t(k1),pairwise.deletion = F))
ggtree(d1 ,ladderize = T)+ geom_point(col="#22AACC")+ theme(legend.position = "right")+geom_tiplab(size=2,col="#22AACC")



E2 = rt[which(rt[,1]=="E2"),]
E2 = E2[,-1]
rownames(E2)=E2[,1]
E2 = E2[,-1]
sort(table(E2$allele1),decreasing = T)
E2_1 = E2[which(E2$allele1=="0"),]
E2_2 = E2[which(E2$allele1=="139:2I"),]
E2_3 = E2[which(E2$allele1=="126:11D"),]
E2_4 = E2[which(E2$allele1=="136:1D"),]
E2_5 = E2[which(E2$allele1=="126:13D"),]
E2_6 = E2[which(E2$allele1=="136:3D"),]


k1=t(E2_1)
d1 = nj(dist.gene(t(k1),pairwise.deletion = F))
ggtree(d1 ,ladderize = T)+ geom_point(col="royalblue")+ theme(legend.position = "right")+geom_tiplab(size=2,col="royalblue")

k1=t(E2_2)
d1 = nj(dist.gene(t(k1),pairwise.deletion = F))
ggtree(d1 ,ladderize = T)+ geom_point(col="#771111")+ theme(legend.position = "right")+geom_tiplab(size=2,col="#771111")

k1=t(E2_3)
d1 = nj(dist.gene(t(k1),pairwise.deletion = F))
ggtree(d1 ,ladderize = T)+ geom_point(col="#22AACC")+ theme(legend.position = "right")+geom_tiplab(size=2,col="#22AACC")

k1=t(E2_4)
d1 = nj(dist.gene(t(k1),pairwise.deletion = F))
ggtree(d1 ,ladderize = T)+ geom_point(col="#666611")+ theme(legend.position = "right")+geom_tiplab(size=2,col="#666611")


E3 = rt[which(rt[,1]=="E3"),]
E3 = E3[,-1]
rownames(E3)=E3[,1]
E3 = E3[,-1]
k3=t(E3)
d3 = nj(dist.gene(t(k3),pairwise.deletion = F))
ggtree(d3)+ geom_point()
ggtree(d3 ,ladderize = T)+ geom_point()+ theme(legend.position = "right")+geom_tiplab(size=2)


E5 = rt[which(rt[,1]=="E5"),]
E5 = E5[,-1]
rownames(E5)=E5[,1]
E5 = E5[,-1]
k5=t(E5)
d5 = nj(dist.gene(t(k5),pairwise.deletion = F))
ggtree(d5 ,ladderize = T)+ geom_point()+ theme(legend.position = "right")+geom_tiplab(size=2)


E6 = rt[which(rt[,1]=="E6"),]
E6 = E6[,-1]
rownames(E6)=E6[,1]
E6 = E6[,-1]
k6=t(E6)
d6 = nj(dist.gene(t(k6),pairwise.deletion = F))
ggtree(d6 ,ladderize = T)+ geom_point()+ theme(legend.position = "right")+geom_tiplab(size=2)




E7 = rt[which(rt[,1]=="E7"),]
E7 = E7[,-1]
rownames(E7)=E7[,1]
E7 = E7[,-1]
k7=t(E7)
d7 = nj(dist.gene(t(k7),pairwise.deletion = F))
ggtree(d7)+ geom_point()
ggtree(d7 ,ladderize = T)+ geom_point()+ theme(legend.position = "right")+geom_tiplab(size=2)

