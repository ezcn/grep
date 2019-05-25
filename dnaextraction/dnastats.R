myd=read.table("/home/enza/oogaProtocol/IMMA/abortion_db/fromValentina/RowDataMiscarriageDNAnames.tsv", header=T , na =NA, sep="\t"  )


mycol=c("#00AFBB", "#E7B800", "#FC4E07")

pDNA<- ggplot(myd, aes(type, (ng.ul*Vol.ul)/1000,  color=Extraction ) ) +geom_boxplot() +facet_grid(Type ~ . , scales="free_y") +scale_color_manual(values=mycol ) +scale_fill_manual (values=mycol ) +theme_bw()  +ylab ("Total Amount of DNA per PoC (ug)" )  +ggtitle("DNA extraction from PoC" )
ggsave("DNAextraction.png", plot= pDNA, device="png", width = 20, height = 15, units = "cm", dpi = 300)

#silvia metter micro invece di u - cambiare induced ...VTP PL - nome legenda "DNA extraction kit"  
