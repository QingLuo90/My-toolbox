library("biomaRt")
listMarts()
##                biomart               version
## 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 98
## 2   ENSEMBL_MART_MOUSE      Mouse strains 98
## 3     ENSEMBL_MART_SNP  Ensembl Variation 98
## 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 98

##The useMart() function can now be used to connect to a specified BioMart database, this must be a valid name given by listMarts(). In the next example we choose to query the Ensembl BioMart database.
ensembl <- useMart("ensembl")
#see what datasets are available
attributes = listAttributes(ensembl)

ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
#or ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

searchAttributes(mart = ensembl)
searchAttributes(mart = ensembl, pattern = "mgi")
#for human gene it hgnc; for mouse gene it is mgi

#read in data
rt = read.table("id.txt",header=T)

id = getBM(attributes = c('ensembl_gene_id', 'mgi_symbol', 'chromosome_name',
                     'start_position', 'end_position', 'band'),
      filters = 'ensembl_gene_id', 
      values = rt[,1], 
      mart = ensembl)
