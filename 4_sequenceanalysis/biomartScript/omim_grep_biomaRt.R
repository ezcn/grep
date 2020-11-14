# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
#####

listMarts()

ensembl <- biomaRt::useMart(host = "http://www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

allGeneOMIM1<-getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","external_gene_name","mim_gene_accession","mim_gene_description"), filters = "", values = "", ensembl, curl = NULL,checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE,quote = "\"")

allGeneOMIM2<-getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","external_gene_name","mim_morbid_description","mim_morbid_accession"), filters = "", values = "", ensembl, curl = NULL,checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE,quote = "\"")

allGeneOMIM_def = merge(allGeneOMIM1,allGeneOMIM2,by=c("chromosome_name","start_position","end_position","ensembl_gene_id","external_gene_name"))

grepOutput<-read.table("/Users/gianlucadamaggio/projects/miscarriage/pathway/grep/20201028_grep_analysis/data/AllGrep.totally_filtered_variants.tsv", header=T, sep="\t")

geneID = grepOutput %>% select(SYMBOL) %>% distinct()

colnames(geneID)[colnames(geneID)=="SYMBOL"] <- "external_gene_name"

grepOMIM=merge(allGeneOMIM_def, geneID, by="external_gene_name")

write.table(grepOMIM, "grepOMIM.tsv", row.names = FALSE, quote = FALSE, col.names = T)
