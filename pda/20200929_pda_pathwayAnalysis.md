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

mydS=read.table("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20200929/data/mergedPda_hgdpIgenomix.tsv", header=T , sep="\t")
mydS$SYMBOL=as.character(mydS$SYMBOL)

geneID = mydS %>% selectd(SYMBOL) %>% distinctd()
hs <- org.Hs.eg.db

#### Transform column in vector:
my.symbols = geneID[,1]
sample=select(hs,keys = as.character(my.symbols), columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
write.table(sample, "/Users/gianlucadamaggio/projects/memorial_exome/pathway/20200929/pda_geneID_to_entrezID.tsv" ,sep="\t", row.names=F)
entrezID=sample$ENTREZID

gene=entrezID
yy = enrichPathway(gene, organism = "human", pvalueCutoff = 0.5, pAdjustMethod = "BH", qvalueCutoff = 0.5, minGSSize = 10, maxGSSize = 500, readable = FALSE)
recap=head(as.data.frame(yy))

write.table(recap, "/Users/gianlucadamaggio/projects/memorial_exome/pathway/20200929/aggregate/recap/pda_aggregate_recap.tsv" ,sep="\t", row.names=F)
write.table(yy, "/Users/gianlucadamaggio/projects/memorial_exome/pathway/20200929/aggregate/pda_reactomeAnalisys.tsv" ,sep="\t", row.names=F)

emapplot(yy, color="pvalue")
ggsave("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20200929/plot/pda_aggregate_emap.png")

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

mydS=read.table("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20200929/data/mergedPda_hgdpIgenomix.tsv", header=T , sep="\t")
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

  write.table(recap, paste("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20200929/individual_analysis/recap/pda_recap_sample_",i,".tsv",sep='') ,sep="\t", row.names=F)
  write.table(yy, paste("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20200929/individual_analysis/pda_allPathways_sample_",i,".tsv",sep=''),sep="\t", row.names=F)

  emapplot(yy, color="pvalue")
  ggsave(paste("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20200929/plot/pda_emap_sample_",i,".png", sep=''))

}

```

### Gene name for each pathway
```
cat pda_aggregate_recap.tsv | cut -f 8 | grep -v gene | tr "/" " " | tr -dc '[:alnum:],[:blank:]\r\n' | tr "\n" ","

for i in 54998 8664 8667 10480 23395 29088 54948 219927 6150 11222 64981 64978 26589 51642 116540 116541 78988 3396 51116 54460 51649 28957 65993 64969 51081 5018 55037 6134 6136 23521 9045 200916 6155 6158 11224 6223 6734 51067 ; do cat ../../pda_geneID_to_entrezID.tsv | grep -w "$i" | cut -f 1 | tr "\n" "," ; done
```
