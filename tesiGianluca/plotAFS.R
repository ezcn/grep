library(ggplot2) 
library(tidyverse)

memory.limit(size=64000)


genome <- read.table("/mpba0/vcolonna/gianluca/TESI/AFS/MergedTSV/TSVforR/merged.chr1.AFS.fb.vep.tsv",header=T,sep="\t")

for (ch in seq(2,22)){
   
   folder="/mpba0/vcolonna/gianluca/TESI/AFS/MergedTSV/TSVforR/"
   
   filename=paste(folder, "merged.chr", ch, ".AFS.fb.vep.tsv", sep="")
   
   myd=read.table(filename, header=T, sep="\t" )
   genome=rbind(genome, myd)
  
}


mysummary <- genome %>% group_by(CSQrank) %>% tally()                                                                                         
write.table(mysummary, file = "/mpba0/vcolonna/gianluca/genome_summaryCSQrank.tsv", append = FALSE, quote = FALSE, sep= "\t") 

myplot = ggplot(subset(genome, CSQfreq>0 , !is.na(CSQfreq)), aes(x =CSQfreq,  fill=Consequence))  + geom_bar(stat = "count", position = "dodge" , width = 0.02) + scale_y_log10() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ ggtitle("AFS") + facet_wrap( ~ CSQrank , ncol=2, scales="free_y" )

ggsave("/mpba0/vcolonna/gianluca/AFSgenome.CSQfreq.dodge.pdf", plot= myplot, width = 50, height = 25, units = "cm", dpi = 300)


genome$CSQrank <- factor(genome$CSQrank, levels = c("MODIFIER","LOW","MODERATE","HIGH"))                       

myplot = ggplot(subset(genome, MAF>0 ), aes(x =MAF,  fill=CSQrank))  + geom_bar(stat = "count" , width = 0.02) + scale_y_log10() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ ggtitle("AFS") + theme_bw() + scale_fill_brewer()   

ggsave("/mpba0/vcolonna/gianluca/AFSgenome.MAF.pdf", plot= myplot, width = 50, height = 25, units = "cm", dpi = 300)

myplot = ggplot(subset(genome, MAF>0 ), aes(x =MAF,  fill=CSQrank))  + geom_bar(stat = "count", position="dodge", width = 0.02) + scale_y_log10() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ ggtitle("AFS") + theme_bw() + scale_fill_brewer()

ggsave("/mpba0/vcolonna/gianluca/AFSgenome.MAF.dodge.pdf", plot= myplot, width = 50, height = 25, units = "cm", dpi = 300)



