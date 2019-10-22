library(ggplot2)                                                                                                                              
library(dplyr)                                                        
library(tidyverse) 
library(cowplot)
library(ggplot2)
library(dplyr)   
library(tidyverse) 

memory.limit(size=32000)

grep<- read.table("/mpba0/vcolonna/gianluca/TESI/hgdp/bothPos/grep.chr22.PyProcessed.tsv" , header=T, sep="\t")
hgdp<- read.table("/mpba0/vcolonna/gianluca/TESI/hgdp/bothPos/hgdpPyProcessed.tsv" , header=T, sep="\t")

hgdp$CSQfreq <- as.numeric(as.character(hgdp$CSQfreq))
grep$CSQfreq <- as.numeric(as.character(grep$CSQfreq))

### uso tutti i cicli
hgdpSel<-select(hgdp,Chr,Pos,VariantClass,CSQallele,CSQrank,Consequence,CSQfreq,REFfreq,ALTfreq,MAF,Population)
full <- rbind(hgdpSel,grep)

### scelgo solo 1 ciclo in particolare
#hgdpSel<-select(hgdp,Chr,Pos,VariantClass,CSQallele,CSQrank,Consequence,CSQfreq,REFfreq,ALTfreq,MAF,Cycle,Population)
#hgdpFiltC1 <- hgdpSel %>% filter ( Cycle == 1 )                                                                                               
#hgdpC1<-select(hgdpFiltC1,Chr,Pos,VariantClass,CSQallele,CSQrank,Consequence,CSQfreq,REFfreq,ALTfreq,MAF,Population)                          
#full <- rbind(hgdpC1, grep) 


########


myplot<-ggplot(subset(full, !is.na(CSQfreq)) , aes(x=Population, y=CSQfreq, fill= VariantClass )) + geom_boxplot(position = position_dodge(preserve = "single")) + facet_wrap( ~ Consequence)


ggsave("/mpba0/vcolonna/gianluca/TESI/hgdp/bothPos/plotWrapConsequence.png", plot = myplot, dpi=300, units="cm", width=50, height =40)

myplot2<-ggplot(subset(full, CSQfreq > 0,  !is.na(CSQfreq)) , aes(x=Population, y=CSQfreq, fill= VariantClass )) + geom_boxplot(position = position_dodge(preserve = "single")) + facet_wrap( ~ Consequence)

ggsave("/mpba0/vcolonna/gianluca/TESI/hgdp/bothPos/plotWrapConsequenceCSQfreqMTZ.png", plot = myplot2, dpi=300, units="cm", width=50, height =40)

myplot3<-ggplot(subset(full, !is.na(CSQfreq)) , aes(x=Population, y=CSQfreq, fill= VariantClass )) + geom_boxplot(position = position_dodge(preserve = "single")) + facet_wrap( ~ CSQrank)

ggsave("/mpba0/vcolonna/gianluca/TESI/hgdp/bothPos/plotWrapCSQrank.png", plot = myplot3, dpi=300, units="cm", width=50, height =40)

myplot4<-ggplot(subset(full,CSQfreq>0, !is.na(CSQfreq)) , aes(x=Population, y=CSQfreq, fill= VariantClass )) + geom_boxplot(position = position_dodge(preserve = "single")) + facet_wrap( ~ CSQrank)

ggsave("/mpba0/vcolonna/gianluca/TESI/hgdp/bothPos/plotWrapCSQrankMTZ.png", plot = myplot4, dpi=300, units="cm", width=50, height =40)
