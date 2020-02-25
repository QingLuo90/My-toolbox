setwd("C:\\Users\\think\\Desktop")
m <- read.table("Taiwan_bases_context_sorted.txt",header=T,stringsAsFactors = F)

x <- cbind(  as.data.frame(table(m$B1))[,2], as.data.frame(table(m$B2))[,2], as.data.frame(table(m$B3))[,2] , as.data.frame(table(m$B4))[,2], as.data.frame(table(m$B5))[,2],as.data.frame(table(m$B6))[,2], as.data.frame(table(m$B7))[,2], as.data.frame(table(m$B8))[,2], as.data.frame(table(m$B9))[,2], as.data.frame(table(m$B10))[,2])

y <- cbind(  as.data.frame(table(m$B12))[,2], as.data.frame(table(m$B13))[,2], as.data.frame(table(m$B14))[,2] , as.data.frame(table(m$B15))[,2], as.data.frame(table(m$B16))[,2],as.data.frame(table(m$B17))[,2], as.data.frame(table(m$B18))[,2], as.data.frame(table(m$B19))[,2], as.data.frame(table(m$B20))[,2], as.data.frame(table(m$B21))[,2])

z <- c(0,0,0,sum(x[,1]))
z <- cbind(x,z,y)
pcm <- z

suppressPackageStartupMessages(library(motifStack))


rownames(pcm) <- c("A","C","G","T")
motif <- new("pcm", mat=as.matrix(pcm), name="MutationContext")
Sys.setenv(R_GSCMD=file.path("C:", "Program Files", "gs", 
                             "gs9.23", "bin", "gswin64c.exe"))

par(mai=c(1,1,1,1))
plot(motif, ic.scale=FALSE, ylab="probability")