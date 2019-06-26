if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)

#####

listMarts()

ensembl <- biomaRt::useMart(host = "http://www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

filters = listFilters(ensembl)                 
attributes = listAttributes(ensembl)


allGene<-getBM(attributes = c("ensembl_gene_id", "chromosome_name","strand","description","start_position","end_position"), filters = "", values = "", ensembl, curl = NULL,checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE,quote = "\"")
##allGene<-filter(prova, ensembl_gene_id=="ENSG00000241186") ## voglio vedere se trova un id a caso.

idGeneEmbryo<-read.table("/home/gianluca/geneidENS.txt")

colnames(idGeneEmbryo)[colnames(idGeneEmbryo)=="V1"] <- "ensembl_gene_id"

startEndEmbryoGene<-merge(idGeneEmbryo,allGene)

write.table(startEndEmbryoGene, "geneset_GO0009790_embryodevelopment.tsv", row.names = FALSE, quote = FALSE)
write.table(allGene, "allGeneStartEndPosition.tsv", row.names = FALSE, quote = FALSE)
