library(maftools)
 setwd( "C:\\Users\\think\\Desktop\\52new")
 maf = read.maf("oncotator_sig22_at.maf")



drive = oncodrive(maf, AACol = NULL, minMut = 3, pvalMethod = "zscore",
                  nBgGenes = 100, bgEstimate = TRUE, ignoreGenes = NULL)
write.table(drive,"new_drive.txt",sep="\t")

colors = c( "#771111", "#117755",'#4477AA',"#771144","#DDDD77","#777711")

lollipopPlot(maf,gene=as.vector(drive$Hugo_Symbol)[1],AACol = "AAChange",pointSize = 3, domainColors = colors,domainLabelSize=3,collapsePosLabel=T,repel=T)

col <- c( "#005594", "#9A2464",'#1EAA83',"#F9D800")

names(col) =c('Splice_Site','Missense_Mutation',"Multi_Hit",'Nonsense_Mutation')

oncoplot(maf = maf, top = 11,genesToIgnore="TTN",legendFontSize = 10,annotationFontSize = 10,annotationTitleFontSize = 10)

oncostrip(maf = maf, genes = as.vector(drive$Hugo_Symbol)[1:15],colors=col)

somaticInteraction = somaticInteractions(maf,top=25)


maf.pfam = pfamDomains(maf, top=10)
maf.pfam$proteinSummary


formutsig = prepareMutSig(maf=maf)
write.table(formutsig,"0.52mutsig.txt",sep="\t",quote = F)




