myd=read.table("/home/enza/oogaProtocol/IMMA/abortion_db/fromValentina/RowDataMiscarriageDNAnames.tsv", header=T , na =NA, sep="\t"  )


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

> myd %>% group_by(Type, type  )  %>% tally() 
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


mydType <- myd %>% group_by( Type, type) %>% summarize(nb=length(type) , max = max (ng.ul*Vol.ul)/1000, min=min(ng.ul*Vol.ul)/1000, average=mean((ng.ul*Vol.ul)/1000), stdev=sd((ng.ul*Vol.ul)/1000)  )

pYield <-  ggplot(mydType, aes(Type, average , color=type   )  ) +geom_point(aes(size=nb )) +geom_errorbar( aes(ymax = average + stdev, ymin=average - stdev, width=0.1) )+scale_color_manual(values=mycol ) +theme_bw() +xlab("" ) + ylab("DNA yield from PoC (ug, average and C.I.)" )  + ggtitle("Tissue homogenization")
ggsave("DNAyield.png", plot= pYield, device="png", width = 15, height = 10, units = "cm", dpi = 300)

#### SILVIA  titololegenda /nb/PoC sample size/  e /type//   ti spiego a voce  


> myd %>% group_by(Extraction  )  %>% tally() 
# A tibble: 3 x 2
  Extraction     n
  <fct>      <int>
1 InstaGene™    52
2 Nucleon       43
3 QIAmp®       146


> myd %>% group_by(Extraction , type  )  %>% tally() 
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


mydExtract <- myd %>% group_by( Extraction, type) %>% summarize(nb=length(type) , max = max (ng.ul*Vol.ul)/1000, min=min(ng.ul*Vol.ul)/1000, average=mean((ng.ul*Vol.ul)/1000), stdev=sd((ng.ul*Vol.ul)/1000)  )
pExtract <- ggplot(mydExtract, aes(Extraction, average , color=type   )  ) +geom_point(aes(size=nb )) +geom_errorbar( aes(ymax = average + stdev, ymin=average - stdev, width=0.1) )+scale_color_manual(values=mycol ) +theme_bw() +xlab("" ) + ylab("DNA yield from PoC (ug, average and C.I.)" )  + ggtitle("Me lo deve dire Vale") 
ggsave("DNAyieldKit.png", plot= pExtract, device="png", width = 15, height = 10, units = "cm", dpi = 300)


#~~~~~~~~~~~~~~~~ make panel 

library(gridExtra)
library (grid) 
library(lattice) 

lay <- rbind(c(1,2)) 
             
myplot<- grid.arrange(pYield, pExtract, nrow = 1, layout_matrix = lay)
ggsave("panelDNA.png", plot = myplot, dpi=300, units="cm", width=30, height =10)


