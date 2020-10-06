library(magrittr)
library(org.Hs.eg.db)
library(ggplot2)
library(ReactomePA)
library(annotate)

distinctd=dplyr::distinct
filterd=dplyr::filter
selectd=dplyr::select

mydS=read.table("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20201006/data/mergedPda_hgdpIgenomix.curation.tsv", header=T , sep="\t")
mydS$SYMBOL=as.character(mydS$SYMBOL)

geneID = mydS %>% selectd(SYMBOL) %>% distinctd()
hs <- org.Hs.eg.db

my.symbols = geneID[,1]
sample=select(hs,keys = as.character(my.symbols), columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
write.table(sample, "/Users/gianlucadamaggio/projects/memorial_exome/pathway/20201006/pda_geneID_to_entrezID.tsv" ,sep="\t", row.names=F)
entrezID=sample$ENTREZID

gene=entrezID
yy = enrichPathway(gene, organism = "human", pvalueCutoff = 0.5, pAdjustMethod = "BH", qvalueCutoff = 0.5, minGSSize = 10, maxGSSize = 500, readable = FALSE)
recap=head(as.data.frame(yy))

write.table(recap, "/Users/gianlucadamaggio/projects/memorial_exome/pathway/20201006/aggregate/recap/pda_aggregate_recap.tsv" ,sep="\t", row.names=F)
write.table(yy, "/Users/gianlucadamaggio/projects/memorial_exome/pathway/20201006/aggregate/pda_reactomeAnalisys.tsv" ,sep="\t", row.names=F)

emapplot(yy, color="pvalue")
ggsave("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20201006/plot/pda_aggregate_emap.png")


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

  write.table(recap, paste("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20201006/individual_analysis/recap/pda_recap_sample_",i,".tsv",sep='') ,sep="\t", row.names=F)
  write.table(yy, paste("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20201006/individual_analysis/pda_allPathways_sample_",i,".tsv",sep=''),sep="\t", row.names=F)

  emapplot(yy, color="pvalue")
  ggsave(paste("/Users/gianlucadamaggio/projects/memorial_exome/pathway/20201006/plot/pda_emap_sample_",i,".png", sep=''))

}
