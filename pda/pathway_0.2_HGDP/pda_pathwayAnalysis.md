# Pathway analysis for PDA samples

### Requirements

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

BiocManager::install("annotate")

BiocManager::install("ReactomePA")

```

### Aggregate analysis

```

library(magrittr)
library(org.Hs.eg.db)
library(ggplot2)
library(ReactomePA)
library(annotate)

distinctd=dplyr::distinct
filterd=dplyr::filter
selectd=dplyr::select

mydS=read.table("/Users/gianlucadamaggio/projects/memorial_exome/pathway/data/pda_full_filtered_noLow.tsv", header=T , sep="\t")
mydS$gene_symbol=as.character(mydS$gene_symbol)

geneID = mydS %>% selectd(gene_symbol) %>% distinctd()
hs <- org.Hs.eg.db

#### Transform column in vector:
my.symbols = geneID[,1]
sample=select(hs,keys = as.character(my.symbols), columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
write.table(sample, "/Users/gianlucadamaggio/projects/memorial_exome/pathway/ReactomePA/pda_geneID_to_entrezID.tsv" ,sep="\t", row.names=F)
entrezID=sample$ENTREZID

gene=entrezID
yy = enrichPathway(gene, organism = "human", pvalueCutoff = 0.5, pAdjustMethod = "BH", qvalueCutoff = 0.5, minGSSize = 10, maxGSSize = 500, readable = FALSE)
recap=head(as.data.frame(yy))

write.table(recap, "/Users/gianlucadamaggio/projects/memorial_exome/pathway/ReactomePA/aggregate/recap/pda_aggregate_recap.tsv" ,sep="\t", row.names=F)
write.table(yy, "/Users/gianlucadamaggio/projects/memorial_exome/pathway/ReactomePA/aggregate/pda_reactomeAnalisys.tsv" ,sep="\t", row.names=F)

emapplot(yy, color="pvalue")
ggsave("/Users/gianlucadamaggio/projects/memorial_exome/pathway/ReactomePA/plot/pda_aggregate_emap.png")

```

### Individual analysis

```

library(magrittr)
library(org.Hs.eg.db)
library(ggplot2)
library(ReactomePA)
library(annotate)


distinctd=dplyr::distinct
filterd=dplyr::filter
selectd=dplyr::select

mydS=read.table("/Users/gianlucadamaggio/projects/memorial_exome/pathway/data/pda_full_filtered_noLow.tsv", header=T , sep="\t")
mydS$gene_symbol=as.character(mydS$gene_symbol)
ID = mydS %>% selectd(sample) %>% distinctd()
hs <- org.Hs.eg.db

for(i in ID[,]){

  nam = paste('sample_',i,sep='')
  temp= mydS %>% selectd(sample, gene_symbol) %>% filterd(sample==i) %>% distinctd(gene_symbol)
  to_process=assign(nam, temp)
  sample=select(hs,keys = as.character(to_process[,]), columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
  entrezID=sample$ENTREZID
  gene=entrezID
  yy = enrichPathway(gene, organism = "human", pvalueCutoff = 0.5, pAdjustMethod = "BH", qvalueCutoff = 0.5, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  recap=head(as.data.frame(yy))

  write.table(recap, paste("/Users/gianlucadamaggio/projects/memorial_exome/pathway/ReactomePA/individual_analysis/recap/pda_recap_sample_",i,".tsv",sep='') ,sep="\t", row.names=F)
  write.table(yy, paste("/Users/gianlucadamaggio/projects/memorial_exome/pathway/ReactomePA/individual_analysis/pda_allPathways_sample_",i,".tsv",sep=''),sep="\t", row.names=F)

  emapplot(yy, color="pvalue")
  ggsave(paste("/Users/gianlucadamaggio/projects/memorial_exome/pathway/ReactomePA/plot/pda_emap_sample_",i,".png", sep=''))

}

```
