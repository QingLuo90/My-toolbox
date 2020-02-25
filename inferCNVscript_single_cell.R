#Change this path to set the working directory
setwd("/Users/yingxinlin/HonoursThesis/SJTU/inferCNV/inferCNVcodesForXianbin")
source("inferCNV.R")

library(gtools)
library(ape) # for the write.tree function
library(Imap)

#Reading the cell types information
cluster<-read.delim("Sample list-6 groups.txt",row.names=1,header=TRUE)
rownames(cluster)<-gsub("_big","",rownames(cluster))
table(cluster[,1])

#Reading the tpt information
tpt<-gsub("_.*","",rownames(cluster))
table(tpt)

#Reading the gene counts data
sjtu_raw= read.delim("timecourse_genecounts.txt",skip=1,row.names=1,header=TRUE)
rownames(sjtu_raw)<-gsub("\\..+","",rownames(sjtu_raw))
sjtu_counts = as.matrix(sjtu_raw[,6:ncol(sjtu_raw)])

#Rename the column names of the gene count matrix
colnames(sjtu_counts)<-gsub("Aligned.sortedByCoord.out.bam","",colnames(sjtu_counts))

#Reorder the column of the gene count matrix according to the cell type information
sjtu_counts<-sjtu_counts[,rownames(cluster)]
#sjtu_counts<-sjtu_counts[,tpt=="E16.5"]
dim(sjtu_counts)



#sjtu_coords_raw = sjtu_raw[,1:2]

#Get the gene length
gene_length<-sjtu_raw[,5]
names(gene_length)<-rownames(sjtu_raw)

#Get the Chromosome of each gene
sjtu_chr = sjtu_raw[,1]
sjtu_chr = gsub(";.+","",sjtu_chr)
table(sjtu_chr)

#Get the location of each gene (start position)
sjtu_pos = sjtu_raw[,2]
sjtu_pos = gsub(";.+","",sjtu_pos)

#Combine the Chromosome and Location information into a matrix
sjtu_coords = cbind(sjtu_chr,sjtu_pos)
rownames(sjtu_coords) <- rownames(sjtu_raw)
sjtu_coord = apply(cbind(sjtu_chr,sjtu_pos),1,paste0,collapse="_")
names(sjtu_coord) <- rownames(sjtu_counts)

#Exclude the ERCC genes
keep = !grepl("ERCC",sjtu_coord)
table(keep)
sjtu_counts = sjtu_counts[keep,]
sjtu_coord = sjtu_coord[keep]

#Order the genes according their chromosome and location information
o = mixedorder(sjtu_coord)
sjtu_counts = sjtu_counts[o,]
sjtu_coord = sjtu_coord[o]


#Transform the count matrix into CPM values
sjtu_CPM = apply(sjtu_counts,2,function(x)1e6*x/sum(x))
sjtu_log2CPM = log2(1+sjtu_CPM)


#Reading the bulk matrix
bulk = as.matrix(read.delim("GSE58827_FPKM.txt",header=TRUE,
                            skip = 1,row.names = 1))
bulkEarly = bulk[,1:12] # 3 reps of the first 4 time points
# need to match the rownames from transcript ids to gene name

library(AnnotationDbi)
library(org.Mm.eg.db)


#Mapping the gene ensembl with gene symbol
mapped = select(org.Mm.eg.db, rownames(sjtu_log2CPM), "SYMBOL", "ENSEMBL")

#Find the genes that are common in single-cell and bulk dataset
idsinBoth = mapped[mapped[,2] %in% rownames(bulkEarly),]
# pull out the first non duplicated for these ids
keep = (!duplicated(idsinBoth[,1]) & 
          !duplicated(idsinBoth[,2]) & 
          (idsinBoth[,2] %in% rownames(bulkEarly)))

idsinBothFilt = idsinBoth[keep,]
dim(idsinBothFilt)



sjtu_log2CPM_matched<-sjtu_log2CPM[idsinBothFilt[,1],]
dim(sjtu_log2CPM_matched)

bulkEarly_matched = bulkEarly[idsinBothFilt[,2],]
bulkEarly_log2CPM_matched = log2(1+bulkEarly_matched)
colnames(bulkEarly_log2CPM_matched) <- paste0("bulk_",colnames(bulkEarly_log2CPM_matched))


#######################################################################
####Generate the plots using bulk reference
#######################################################################
data = cbind(sjtu_log2CPM_matched,bulkEarly_log2CPM_matched)

gene_order = matrix(0,ncol=1,nrow=nrow(data))
rownames(gene_order) = rownames(data)
cutoff = 4.5
reference_obs = colnames(data)[grepl("bulk",colnames(data))]
transform_data = FALSE
window_length = 101
max_centered_threshold = 20
noise_threshold = 0.3
num_ref_groups = 1
##Setting the output path
out_path = "/Users/yingxinlin/HonoursThesis/SJTU/inferCNV/inferCNVcodesForXianbin"
plot_steps = FALSE

data_filt = data
gene_order_filt = gene_order

dim(data_filt)

result = infer_cnv(data_filt,
                   gene_order_filt,
                   cutoff,
                   reference_obs,
                   transform_data,
                   window_length,
                   max_centered_threshold,
                   noise_threshold,
                   num_ref_groups,
                   out_path,
                   plot_steps)
# save(result,file = "/Users/yingxinlin/HonoursThesis/SJTU/inferCNV/P3.25/bulkRef/resultwithbulkreference.RData")


sjtu_coords_filt = sjtu_coords

plot_data = result$VIZ

dim(plot_data)

# clip off the first and last 15 rows since the density estimations gets weird at the ends
plot_data_filt <- plot_data[10:(nrow(plot_data)-10),] 

rownames(sjtu_coords_filt) <- rownames(sjtu_coords_filt)
rownames(plot_data_filt)<- rownames(plot_data_filt)


contigs = sjtu_coords_filt[rownames(plot_data_filt),1]
reference_idx = result$REF_OBS_IDX
ref_groups = result$REF_GROUPS
##Setting the output path
out_dir = 
title = "bulk reference"
obs_title = ""
ref_title = ""
ref_contig = NULL

#Using cluster names as labels
label<-cluster[colnames(sjtu_CPM),1]
table(label)

#Using time points as labels
#label<-gsub("_.*","",colnames(sjtu_CPM))
#table(label)

plot_cnv(plot_data_filt,
         contigs,
         reference_idx,
         ref_contig,
         ref_groups,
         out_dir,
         title,
         obs_title,
         ref_title,
         contig_cex=2,
         k_obs_groups=1,
         color_safe_pal=FALSE,
         celllabels=label)






#######################################################################
####Generate the plots using single-cell data itself as reference
#######################################################################
data = sjtu_log2CPM

gene_order = matrix(0,ncol=1,nrow=nrow(data))
rownames(gene_order) = rownames(data)
cutoff = 4.5
reference_obs = colnames(data)
transform_data = FALSE
window_length = 101
max_centered_threshold = 20
noise_threshold = 0.3
num_ref_groups = 1
##Setting the output path
out_path = "/Users/yingxinlin/HonoursThesis/SJTU/inferCNV/inferCNVcodesForXianbin"
plot_steps = FALSE

data_filt = data
gene_order_filt = gene_order

dim(data_filt)

result_sc = infer_cnv(data_filt,
                   gene_order_filt,
                   cutoff,
                   reference_obs,
                   transform_data,
                   window_length,
                   max_centered_threshold,
                   noise_threshold,
                   num_ref_groups,
                   out_path,
                   plot_steps)

sjtu_coords_filt = sjtu_coords

plot_data = result_sc$VIZ

dim(plot_data)

# clip off the first and last 15 rows since the density estimations gets weird at the ends
plot_data_filt <- plot_data[10:(nrow(plot_data)-10),] 

rownames(sjtu_coords_filt) <- rownames(sjtu_coords_filt)
rownames(plot_data_filt)<- rownames(plot_data_filt)


contigs = sjtu_coords_filt[rownames(plot_data_filt),1]
reference_idx = result$REF_OBS_IDX
ref_groups = result$REF_GROUPS
##Setting the output path
out_dir = "/Users/yingxinlin/HonoursThesis/SJTU/inferCNV/inferCNVcodesForXianbin"
title = "single-cell reference"
obs_title = ""
ref_title = ""
ref_contig = NULL

#Using cluster names as labels
label<-cluster[colnames(sjtu_CPM),1]
table(label)

#Using time points as labels
#label<-gsub("_.*","",colnames(sjtu_CPM))
#table(label)


plot_cnv(plot_data_filt,
         contigs,
         reference_idx,
         ref_contig,
         ref_groups,
         out_dir,
         title,
         obs_title,
         ref_title,
         contig_cex=2,
         k_obs_groups=1,
         color_safe_pal=FALSE,
         celllabels=label)




pdf(file="/Users/yingxinlin/HonoursThesis/SJTU/inferCNV/inferCNVcodesForXianbin/cellTypeColour.pdf",height = 2,width = 10)
barplot(rep(1,6),col=rainbow(6),names=c("Endothelial Cell","Erythrocyte","Hepatoblast","Macrophage","Megakaryocyte","Mesenchymal Cell"))
dev.off()

pdf(file="/Users/yingxinlin/HonoursThesis/SJTU/inferCNV/inferCNVcodesForXianbin/timeColour.pdf",height = 2,width = 6)
barplot(rep(1,5),col=rainbow(5),names=c("E11.5","E12.5","E13.5","E14.5","E16.5"))
dev.off()