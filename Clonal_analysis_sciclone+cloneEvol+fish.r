setwd("D:\\0530\\SciClone\\M2")

library(sciClone)
v1 = read.table("M2_T2_vaf.txt", header=F)
cn1 = read.table("M2_T2_CNV.txt",header=F)
v2 = read.table("M2_T4_vaf.txt", header=F)
cn2 = read.table("M2_T4_CNV.txt",header=F)
names = c("M2_T2","M2_T4")
sc = sciClone(vafs=list(v1,v2),             copyNumberCalls=list(cn1,cn2),cnCallsAreLog2=T,   copyNumberMargins=0.35,minimumDepth=50,sampleNames=names[1:2],useSexChrs=FALSE)
writeClusterTable(sc, "M2_T2T4_clusters")
writeClusterSummaryTable(sc, "M2_T2T4_clusters.summary")

sc.plot2d(sc, "M2_T2T4_figure.pdf")





library(clonevol)
x=read.table("cluster.txt",sep="\t",header=T)

clone.colors <- c( "#771111", "#117755",'#4477AA',"#771144","#DDDD77","#777711")

vaf.col.names <- grep('.vaf', colnames(x), value=T)
sample.names <- gsub('.vaf', '', vaf.col.names)
x[, sample.names] <- x[, vaf.col.names]
vaf.col.names <- sample.names
# prepare sample grouping
sample.groups <- c('M2_T2', 'M4_T4');
names(sample.groups) <- vaf.col.names
# setup the order of clusters to display in various plots (later)
x <- x[order(x$cluster),]


pp <- plot.variant.clusters(x,
                            cluster.col.name = 'cluster',
                            show.cluster.size = FALSE,
                            cluster.size.text.color = 'blue',
                            vaf.col.names = vaf.col.names,
                            vaf.limits = 70,
                            sample.title.size = 20,
                            violin = FALSE,
                            box = FALSE,
                            jitter = TRUE,
                            jitter.shape = 16,
                            jitter.color = clone.colors,
                            jitter.size = 1.5,
                            jitter.alpha = .5,
                            jitter.center.method = 'median',
                            jitter.center.size = 1,
                            jitter.center.color = 'darkgray',
                            jitter.center.display.value = 'none',
                            highlight = 'is.driver',
                            highlight.shape = 21,
                            highlight.color = 'blue',
                            highlight.fill.color = 'green',
                            highlight.note.col.name = 'gene',
                            highlight.note.size = 2,
                            order.by.total.vaf = FALSE)


plot.pairwise(x, col.names = vaf.col.names,
              out.prefix = 'variants.pairwise.plot',
              colors = clone.colors)

y = infer.clonal.models(variants = x,
                        cluster.col.name = 'cluster',
                        vaf.col.names = vaf.col.names,
                        sample.groups = sample.groups,
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = 1,
                        cluster.center = 'mean',
                        ignore.clusters = NULL,
                        clone.colors = clone.colors,
                        min.cluster.vaf = 0.01,
                        # min probability that CCF(clone) is non-negative
                        sum.p = 0.05,
                        # alpha level in confidence interval estimate for CCF(clone)
                        alpha = 0.05)

z <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

plot.clonal.models(z,
                   # box plot parameters
                   box.plot = T,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 70,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 25,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 1.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 0.5,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 2.5,
                   mtcab.node.label.size = 0.5,
                   mtcab.node.text.size = .5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = 'output',
                   out.format = 'pdf',
                   overwrite.output = TRUE,width = 9,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,3,3))

dev.off()



library(fishplot)
timepoints=c(0,10,30,250,300)

frac.table = matrix(
    c(5,0,0,0,7,0,0,0,6,3,0,0,44,36,18,4,44,36,18,4),
    ncol=length(timepoints))

parents = c(0,1,2,2)

fish = createFishObject(frac.table,parents,timepoints=timepoints, clone.labels=c("Founding", "Subclone 1","Subclone 2","Subclone 3" ))
fish = layoutClones(fish)
sample.times = timepoints

clone.colors <- c( "#771111", "#117755",'#4477AA',"#771144")
fish = createFishObject(frac.table,parents,timepoints=timepoints, clone.labels=c("Founding", "Subclone 1","Subclone 2","Subclone 3" ),col=clone.colors)
fish = layoutClones(fish)
sample.times = timepoints


fishPlot(fish,shape="spline",title.btm="633734",
         vlines=sample.times, vlab=sample.times, cex.title=0.5,bg.col=c("#771111", "#117755",'#4477AA'))

drawLegend(fish)



library(fishplot)
timepoints=c(0,10,20,200,300)

frac.table = matrix(
    c(5,0,0,0,6,0,0,0,7,2,0,0,44,36,18,4,44,36,18,4),
    ncol=length(timepoints))

parents = c(0,1,2,2)

fish = createFishObject(frac.table,parents,timepoints=timepoints, clone.labels=c("Founding", "Subclone 1","Subclone 2","Subclone 3" ))
fish = layoutClones(fish)
sample.times = timepoints

clone.colors <- c( "#606B74", '#A9A6D3',"#3E1D63","#2C55A7")
fish = createFishObject(frac.table,parents,timepoints=timepoints, clone.labels=c("Founding", "Subclone 1","Subclone 2","Subclone 3" ),col=clone.colors)
fish = layoutClones(fish)
sample.times = timepoints


fishPlot(fish,shape="spline",title.btm="",
         vlines=sample.times, vlab=sample.times, cex.title=0.5,bg.col=c("#117755", "#4477AA",'#771111'),pad.left = 0.05,col.border="black")

drawLegend(fish)
