library(magrittr)
library(org.Hs.eg.db)
library(ggplot2)
library(ReactomePA)
library(annotate)

distinctd=dplyr::distinct
filterd=dplyr::filter
selectd=dplyr::select

mydS=read.table("/Users/gianlucadamaggio/projects/miscarriage/pathway/grep/20201028_PL_analysis/data/PL_totally_filtered_variants.tsv", header=T , sep="\t")
mydS$SYMBOL=as.character(mydS$SYMBOL)

geneID = mydS %>% selectd(SYMBOL) %>% distinctd()
hs <- org.Hs.eg.db

#### Transform column in vector:
my.symbols = geneID[,1]
sample=select(hs,keys = as.character(my.symbols), columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
write.table(sample, "~/projects/miscarriage/pathway/grep/20201028_PL_analysis/PL_geneID_to_entrezID.tsv" ,sep="\t", row.names=F)
entrezID=sample$ENTREZID

gene=entrezID
yy = enrichPathway(gene, organism = "human", pvalueCutoff = 0.06, pAdjustMethod = "BH", qvalueCutoff = 0.06, minGSSize = 10, maxGSSize = 500, readable = TRUE)
recap=head(as.data.frame(yy))

write.table(recap, "~/projects/miscarriage/pathway/grep/20201028_PL_analysis/aggregate/recap/PL_aggregate_recap.tsv" ,sep="\t", row.names=F)
write.table(yy, "~/projects/miscarriage/pathway/grep/20201028_PL_analysis/aggregate/PL_reactomeAnalisys.tsv" ,sep="\t", row.names=F)

emapplot(yy, color="p.adjust")
ggsave("~/projects/miscarriage/pathway/grep/20201028_PL_analysis/plot/PL_aggregate_emap.png")

heatplot(yy)
ggsave("~/projects/miscarriage/pathway/grep/20201028_PL_analysis/plot/PL_aggregate_heat.png")

ID = mydS %>% selectd(ID) %>% distinctd()

for(i in ID[,]){

  nam = paste('sample_',i,sep='')
  temp= mydS %>% selectd(ID, SYMBOL) %>% filterd(ID==i) %>% distinctd(SYMBOL)
  to_process=assign(nam, temp)
  sample=select(hs,keys = as.character(to_process[,]), columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
  entrezID=sample$ENTREZID
  gene=entrezID
  yy = enrichPathway(gene, organism = "human", pvalueCutoff = 0.06, pAdjustMethod = "BH", qvalueCutoff = 0.06, minGSSize = 10, maxGSSize = 500, readable = TRUE)
  recap=head(as.data.frame(yy))

  write.table(recap, paste("~/projects/miscarriage/pathway/grep/20201028_PL_analysis/individual_analysis/recap/PL_recap_sample_",i,".tsv",sep='') ,sep="\t", row.names=F)
  write.table(yy, paste("~/projects/miscarriage/pathway/grep/20201028_PL_analysis/individual_analysis/PL_allPathways_sample_",i,".tsv",sep=''),sep="\t", row.names=F)

  emapplot(yy, color="p.adjust")
  ggsave(paste("~/projects/miscarriage/pathway/grep/20201028_PL_analysis/plot/PL_emap_sample_",i,".png", sep=''))

}
