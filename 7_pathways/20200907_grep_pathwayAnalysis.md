# Pathway analysis for GREP samples

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

mydS=read.table("/Users/gianlucadamaggio/projects/miscarriage/pathway/grep/20200907/data/allGrep_filtered.tsv", header=T , sep="\t")
mydS$SYMBOL=as.character(mydS$SYMBOL)

geneID = mydS %>% selectd(SYMBOL) %>% distinctd()
hs <- org.Hs.eg.db

#### Transform column in vector:
my.symbols = geneID[,1]
sample=select(hs,keys = as.character(my.symbols), columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
write.table(sample, "~/projects/miscarriage/pathway/grep/20200907/grep_geneID_to_entrezID.tsv" ,sep="\t", row.names=F)
entrezID=sample$ENTREZID

gene=entrezID
yy = enrichPathway(gene, organism = "human", pvalueCutoff = 0.5, pAdjustMethod = "BH", qvalueCutoff = 0.5, minGSSize = 10, maxGSSize = 500, readable = FALSE)
recap=head(as.data.frame(yy))

write.table(recap, "~/projects/miscarriage/pathway/grep/20200907/aggregate/recap/grep_aggregate_recap.tsv" ,sep="\t", row.names=F)
write.table(yy, "~/projects/miscarriage/pathway/grep/20200907/aggregate/grep_reactomeAnalisys.tsv" ,sep="\t", row.names=F)

emapplot(yy, color="pvalue")
ggsave("~/projects/miscarriage/pathway/grep/20200907/plot/grep_aggregate_emap.png")

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

mydS=read.table("/Users/gianlucadamaggio/projects/miscarriage/pathway/grep/20200907/data/allGrep_filtered.tsv", header=T , sep="\t")
mydS$SYMBOL=as.character(mydS$SYMBOL)
ID = mydS %>% selectd(ID) %>% distinctd()
hs <- org.Hs.eg.db

for(i in ID[,]){

  nam = paste('sample_',i,sep='')
  temp= mydS %>% selectd(ID, SYMBOL) %>% filterd(ID==i) %>% distinctd(SYMBOL)
  to_process=assign(nam, temp)
  sample=select(hs,keys = as.character(to_process[,]), columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
  entrezID=sample$ENTREZID
  gene=entrezID
  yy = enrichPathway(gene, organism = "human", pvalueCutoff = 0.5, pAdjustMethod = "BH", qvalueCutoff = 0.5, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  recap=head(as.data.frame(yy))

  write.table(recap, paste("~/projects/miscarriage/pathway/grep/20200907/individual_analysis/recap/grep_recap_sample_",i,".tsv",sep='') ,sep="\t", row.names=F)
  write.table(yy, paste("~/projects/miscarriage/pathway/grep/20200907/individual_analysis/grep_allPathways_sample_",i,".tsv",sep=''),sep="\t", row.names=F)

  emapplot(yy, color="pvalue")
  ggsave(paste("~/projects/miscarriage/pathway/grep/20200907/plot/grep_emap_sample_",i,".png", sep=''))

}

```
