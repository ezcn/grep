if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
#####

listMarts()

ensembl <- biomaRt::useMart(host = "http://www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

filters = listFilters(ensembl)                 
attributes = listAttributes(ensembl)

filtersimp<- filter(attributes, name=="imprinted")

allGenePosition<-getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","external_gene_name"), filters = "", values = "", ensembl, curl = NULL,checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE,quote = "\"")
#allGene<-filter(prova, ensembl_gene_id=="ENSG00000241186") ## voglio vedere se trova un id a caso.

idGeneEmbryo<-read.table("/home/gianluca/geneidENS.txt")

colnames(idGeneEmbryo)[colnames(idGeneEmbryo)=="V1"] <- "ensembl_gene_id"

startEndEmbryoGene<-merge(idGeneEmbryo,allGenePosition)
startEndEmbryoGene<-startEndEmbryoGene[c("chromosome_name","start_position","end_position","ensembl_gene_id","external_gene_name")]

write.table(startEndEmbryoGene, "embryodev.bed", row.names = FALSE, quote = FALSE, col.names = F)
write.table(allGene, "allGeneStartEndPosition.tsv", row.names = FALSE, quote = FALSE)
