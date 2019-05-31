myd=read.table("/home/enza/oogaProtocol/IMMA/abortion_db/fromValentina/RowDataMiscarriageDNAnames.tsv", header=T , na =NA, sep="\t"  )

myd$method <-(ifelse(myd$Extraction=="QIAmp®", "Membrane", ifelse(myd$Extraction=="InstaGene™", "Chelex resin",  "Nucleon resin" )))
myd$method<- factor(myd$method, levels=c("Chelex resin", "Nucleon resin", "Membrane"))

myd$tpf <-(ifelse(myd$type=="induced", "VTP", ifelse(myd$type=="miscarriage_first", "FPL",  "RPL" )))
myd$tpf<- factor(myd$tpf, levels=c("VTP", "FPL", "RPL"))

mycol=c("#00AFBB", "#E7B800", "#FC4E07")

pDNA<- ggplot(myd, aes(type, (ng.ul*Vol.ul)/1000,  color=Extraction ) ) +geom_boxplot() +facet_grid(Type ~ . , scales="free_y") +scale_color_manual(values=mycol ) +scale_fill_manual (values=mycol ) +theme_bw()  +ylab ("Total Amount of DNA per PoC (ug)" )  +ggtitle("DNA extraction from PoC" )
ggsave("DNAextraction.png", plot= pDNA, device="png", width = 20, height = 15, units = "cm", dpi = 300)

#silvia metter micro invece di u - cambiare induced ...VTP PL - nome legenda "DNA extraction kit"  

> myd %>% group_by(Type  )  %>% tally() 
# A tibble: 3 x 2
  Type             n
  <fct>        <int>
1 Culture         20
2 Dry            211
3 Lysis Buffer    10
> 

> myd %>% group_by(Type, tpf  )  %>% tally() 
# A tibble: 8 x 3
# Groups:   Type [3]
  Type         type                      n
  <fct>        <fct>                 <int>
1 Culture      miscarriage_first        13
2 Culture      miscarriage_recurrent     7
3 Dry          induced                 127
4 Dry          miscarriage_first        54
5 Dry          miscarriage_recurrent    30
6 Lysis Buffer induced                   6
7 Lysis Buffer miscarriage_first         3
8 Lysis Buffer miscarriage_recurrent     1


#mydType <- subset(myd, Type!="Culture") %>% group_by( Type, type) %>% summarize(nb=length(type) , max = max (ng.ul*Vol.ul)/1000, min=min(ng.ul*Vol.ul)/1000, average=mean((ng.ul*Vol.ul)/1000), stdev=sd((ng.ul*Vol.ul)/1000)  )

###silvia 
mydType <- subset(myd, Type!="culture") %>% group_by( Type, tpf) %>% summarize(nb=length(tpf) , max = max (ng.ul*Vol.ul)/1000, min=min(ng.ul*Vol.ul)/1000, average=mean((ng.ul*Vol.ul)/1000), stdev=sd((ng.ul*Vol.ul)/1000)  )

#pYield <-  ggplot(mydType, aes(Type, average , color=Type   )  ) +geom_point(aes(size=nb )) +geom_errorbar( aes(ymax = average + stdev, ymin=average - stdev, width=0.1) )+scale_color_manual(values=mycol ) +theme_bw() +xlab("" ) + ylab("DNA yield from PoC (ug)" )  + ggtitle("Tissue homogenization") + labs(color= "") + scale_size(name="PoC sample size")  
#ggsave("DNAyield.png", plot= pYield, device="png", width = 15, height = 10, units = "cm", dpi = 300)
###silvia 
pYield<-ggplot(mydType, aes(Type, average , color=tpf   )  ) +geom_point(aes(size=nb )) +geom_errorbar( aes(ymax = average + stdev, ymin=average - stdev, width=0.1) )+scale_color_manual(values=mycol ) +theme_bw() +xlab("" ) + ylab("DNA yield from PoC (microgram)" )  + ggtitle("Tissue homogenization") + labs(color= "") + scale_size(name="PoC sample size")
ggsave("DNAyield.png", plot= pYield, device="png", width = 15, height = 10, units = "cm", dpi = 300)

#, average and C.I.
#### SILVIA  titololegenda /nb/PoC sample size/  e /type//   ti spiego a voce  


> myd %>% group_by(Extraction  )  %>% tally() 
# A tibble: 3 x 2
  Extraction     n
  <fct>      <int>
1 InstaGene™    52
2 Nucleon       43
3 QIAmp®       146


> myd %>% group_by(method , tpf  )  %>% tally() 
# A tibble: 9 x 3
# Groups:   Extraction [3]
  Extraction type                      n
  <fct>      <fct>                 <int>
1 InstaGene™ induced                  24
2 InstaGene™ miscarriage_first        14
3 InstaGene™ miscarriage_recurrent    14
4 Nucleon    induced                  22
5 Nucleon    miscarriage_first        14
6 Nucleon    miscarriage_recurrent     7
7 QIAmp®     induced                  87
8 QIAmp®     miscarriage_first        42
9 QIAmp®     miscarriage_recurrent    17


#mydExtract <- myd %>% group_by( Extraction, type) %>% summarize(nb=length(type) , max = max (ng.ul*Vol.ul)/1000, min=min(ng.ul*Vol.ul)/1000, average=mean((ng.ul*Vol.ul)/1000), stdev=sd((ng.ul*Vol.ul)/1000)  )
###silvia 
mydExtract <- myd %>% group_by( method, tpf) %>% summarize(nb=length(tpf) , max = max (ng.ul*Vol.ul)/1000, min=min(ng.ul*Vol.ul)/1000, average=mean((ng.ul*Vol.ul)/1000), stdev=sd((ng.ul*Vol.ul)/1000)  )

#pExtract <- ggplot(mydExtract, aes(Extraction, average , color=type   )  ) +geom_point(aes(size=nb )) +geom_errorbar( aes(ymax = average + stdev, ymin=average - stdev, width=0.1) )+scale_color_manual(values=mycol ) +theme_bw() +xlab("" ) + ylab("DNA yield from PoC (ug) " )  + ggtitle("DNA isolation") 
#ggsave("DNAyieldKit.png", plot= pExtract, device="png", width = 15, height = 10, units = "cm", dpi = 300)
#, average and C.I.
##silvia 
pExtract <- ggplot(mydExtract, aes(method, average , color=tpf   )  ) +geom_point(aes(size=nb )) +geom_errorbar( aes(ymax = average + stdev, ymin=average - stdev, width=0.1) )+scale_color_manual(values=mycol ) +theme_bw() +xlab("" ) + ylab("DNA yield from PoC (microgram) " )  + ggtitle("DNA isolation") + labs(color= "") + scale_size(name="PoC sample size")
ggsave("DNAyieldKit.png", plot= pExtract, device="png", width = 15, height = 10, units = "cm", dpi = 300)
#~~~~~~~~~~~~~~~~ make panel 

#library(gridExtra)
#ibrary (grid) 
#library(lattice) 

#lay <- rbind(c(1,2)) 
             
#myplot<- grid.arrange(pYield, pExtract, nrow = 1, layout_matrix = lay)
#ggsave("panelDNA.png", plot = myplot, dpi=300, units="cm", width=25, height =10)

##silvia
library(ggpubr)
myplot<-ggarrange(pYield,pExtract, labels = c("A", "B"), common.legend = TRUE, legend = "right")
ggsave("panelDNA.png", plot = myplot, dpi=300, units="cm", width=25, height =10)

