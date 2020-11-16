# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("biomaRt")

# Load Library
library(biomaRt)
library(dplyr)

# BioMart
listMarts()

ensembl <- biomaRt::useMart(host = "http://www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

# Only 2 External attributes can be queried at same time:
# Get "mim_gene_accession","mim_gene_description" from all genes
allGeneOMIM1<-getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","external_gene_name","mim_gene_accession","mim_gene_description"), filters = "", values = "", ensembl, curl = NULL,checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE,quote = "\"")

# Get "mim_morbid_description","mim_morbid_accession" from all genes
allGeneOMIM2<-getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","external_gene_name","mim_morbid_accession","mim_morbid_description"), filters = "", values = "", ensembl, curl = NULL,checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE,quote = "\"")

# Merge two subset for obtain all 4 OMIM External attributes
allGeneOMIM_def = merge(allGeneOMIM1,allGeneOMIM2,by=c("chromosome_name","start_position","end_position","ensembl_gene_id","external_gene_name"))

# Load GREP output Table
grepOutput<-read.table("/Users/gianlucadamaggio/projects/miscarriage/data/grep_output/Grep_final_allfilters.tsv", header=T, sep="\t")

# Select only uniq Gene Name
geneID = grepOutput %>% select(SYMBOL) %>% distinct()

# Change column name for match with BioMart Query
colnames(geneID)[colnames(geneID)=="SYMBOL"] <- "external_gene_name"

# Merge BioMart Query with GREP's output table
grepOMIM=merge(allGeneOMIM_def, geneID, by="external_gene_name")

# Arrange column
grepOMIM= grepOMIM[c(1,2,3,4,5,6,7,9,8)]

# Write table with omim info for each Gene Name of GREP's output
write.table(grepOMIM, "grepOMIM.tsv", row.names = FALSE, quote = F, col.names = T, sep="\t")
